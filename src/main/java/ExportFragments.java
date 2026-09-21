import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scan.StorageStrategy;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.mzml.MZMLFile;
import umich.ms.glyco.Glycan;
import umich.ms.glyco.GlycanParser;
import umich.ms.glyco.GlycanResidue;

import java.io.*;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicBoolean;

public class ExportFragments {

    private int threadNum = 1;
    private ResultProcessor resultProcessor;

    /** Backbone ion types named on the command line; empty means "use the search's own". */
    private final List<String> overrideIonTypes;
    private final boolean addNeutralLoss;

    /** The search's own parameters, read from fragpipe.workflow. */
    private SearchParams search;
    /** Backbone series to generate: the search's standard letters plus every custom definition. */
    private List<FragmentAnnotator.Series> backboneSeries = new ArrayList<>();
    /** Series that carry fragment remainders (msfragger.labile_fragment_ion_series). */
    private List<FragmentAnnotator.Series> labileSeries = new ArrayList<>();
    /** The ion categories this result gets columns for. */
    private List<IonCategory> categories = new ArrayList<>();

    private final String glycanResiduesPath;
    private final String glycanModsPath;
    private HashMap<String, GlycanResidue> glycanResiduesMap = new HashMap<>();

    /** Warn at most once per run that a glyco search has no glycan compositions to name ions with. */
    private final AtomicBoolean warnedNoGlycanColumn = new AtomicBoolean(false);

    public ExportFragments(File resultsFolder, int threadNum, List<String> overrideIonTypes,
            boolean addNeutralLoss, String glycanResiduesPath, String glycanModsPath) throws IOException {
        this.threadNum = threadNum;
        this.overrideIonTypes = overrideIonTypes;
        this.addNeutralLoss = addNeutralLoss;
        this.glycanResiduesPath = glycanResiduesPath;
        this.glycanModsPath     = glycanModsPath;

        importData(resultsFolder);
    }

    private void importData(File resultsFolder) throws IOException {

        resultProcessor = new ResultProcessor(resultsFolder);

        if (resultProcessor.manifestFile != null) {

            if (resultProcessor.psmIndexToName.containsValue("ionint") && resultProcessor.psmIndexToName.containsValue("ionmz")) {
                System.err.println("This psm.tsv already carries ion annotation columns "
                        + "(ions / ion_mz / ion_int), so it was written by a previous run of this tool. "
                        + "Annotation rewrites psm.tsv in place and cannot be repeated on its own output. "
                        + "Re-run Philosopher to regenerate psm.tsv, then annotate it again.");
                System.exit(1);
            }

            configureFromSearch(resultsFolder);

            // Initialise the glycan residue/mod database before spawning threads so that
            // glycoShortNames is fully populated before any concurrent access.
            initGlycanResiduesMap();

            try {
                ExecutorService executorService = Executors.newFixedThreadPool(threadNum);

                for (String expNum : resultProcessor.resultsDict.keySet()) {
                    ArrayList<String[]> onePSMData = readOnePSM(expNum);
                    if (onePSMData == null || onePSMData.isEmpty()) continue;

                    HashMap<String, ArrayList<Integer>> fileToIndices = groupPSMsByFile(onePSMData);

                    LinkedHashMap<IonCategory, ArrayList<IonMatch>[]> ionMatchesMap = buildIonMatchesMap(onePSMData.size());
                    LinkedHashMap<IonCategory, ArrayList<IonMatch>[]> pairedIonMatchesMap =
                            resultProcessor.hasPairedScanNum ? buildIonMatchesMap(onePSMData.size()) : null;

                    // Open all mzML files and queue per-PSM tasks without waiting between files.
                    // Workers fetch spectra concurrently (overlapping I/O with CPU work).
                    ArrayList<Future<?>> futures = new ArrayList<>();
                    for (String fileName : fileToIndices.keySet()) {
                        ScanCollectionDefault scans = openMzmlFile(fileName);
                        if (scans == null) continue;
                        ArrayList<Integer> indices = fileToIndices.get(fileName);
                        System.out.println("Annotating " + indices.size() + " PSMs from " + fileName);
                        for (int psmIndex : indices) {
                            futures.add(executorService.submit(
                                    getOneAnnotation(psmIndex, onePSMData,
                                            ionMatchesMap, pairedIonMatchesMap, scans)));
                        }
                    }
                    for (Future<?> future : futures) future.get();

                    writeOneExperiment(expNum, onePSMData, ionMatchesMap, pairedIonMatchesMap);
                }

                executorService.shutdown();
            } catch (IOException | InterruptedException | ExecutionException e) {
                System.exit(1);
                throw new RuntimeException(e);
            }
            System.exit(0);
        }
    }

    /**
     * Read the search's own settings and decide what to annotate with them: which backbone series,
     * at what tolerance, and whether the three labile categories apply at all.
     */
    private void configureFromSearch(File resultsFolder) {
        search = SearchParams.read(resultsFolder);
        FragmentAnnotator.configureTolerance(search.fragTol, search.fragTolDa);
        System.out.println("Fragment tolerance: " + search.fragTol + (search.fragTolDa ? " Da" : " ppm"));

        // Standard letters come from the search unless the command line named its own. Custom
        // series are always generated: MSFragger generates every definition it is given, whether
        // or not fragment_ion_series also names it.
        List<String> letters = overrideIonTypes.isEmpty()
                ? standardLetters(search.ionSeries)
                : overrideIonTypes;
        if (letters.isEmpty()) {
            System.err.println("WARNING: the search declares no standard ion series; defaulting to b, y");
            letters = Arrays.asList("b", "y");
        }
        List<String> backboneLabels = new ArrayList<>(letters);
        for (CustomIon ci : search.customIons) {
            if (!backboneLabels.contains(ci.name)) backboneLabels.add(ci.name);
        }
        backboneSeries = FragmentAnnotator.resolveSeries(backboneLabels, search.customIons);
        System.out.println("Backbone ion series: " + backboneLabels);

        if (search.labile) {
            // The base series for fragment remainders are the ones MSFragger generated them for,
            // intersected with what is being annotated. Only a search that declared NO labile
            // series falls back to all of them: an empty intersection means the declared series
            // are simply not being annotated, and widening it there would put remainders on series
            // the search never generated any for.
            List<String> labileLabels = new ArrayList<>();
            if (search.labileSeries.isEmpty()) {
                labileLabels = backboneLabels;
            } else {
                for (String l : search.labileSeries) if (backboneLabels.contains(l)) labileLabels.add(l);
                if (labileLabels.isEmpty()) {
                    System.err.println("WARNING: the search generated fragment remainders for "
                            + search.labileSeries + ", none of which are being annotated, so the "
                            + "fragment remainder column will be empty.");
                }
            }
            labileSeries = FragmentAnnotator.resolveSeries(labileLabels, search.customIons);

            if (!search.offsets.isEmpty()) {
                System.out.println("Labile mode: " + search.offsets.size() + " mass offset(s), "
                        + "remainder series " + labileLabels);
            } else if (search.globalOffset != null) {
                // No per-mass entry to match, but the search's global ion lists say what a labile
                // modification produces. This is the shape of a glyco search run without detailed
                // offsets, where the glycan masses live only in the detailed parameter.
                System.out.println("Labile mode: no per-mass offsets; using the search's global "
                        + "diagnostic / Y-type / remainder lists, remainder series " + labileLabels);
            } else {
                System.err.println("WARNING: labile search mode is on but no labile mass offsets "
                        + "could be read from fragpipe.workflow; the labile columns will be empty.");
            }
        }

        categories.add(IonCategory.BACKBONE);
        if (search.labile) {
            categories.add(IonCategory.FRAG_REMAINDER);
            categories.add(IonCategory.PEP_REMAINDER);
            categories.add(IonCategory.DIAGNOSTIC);
        }
    }

    /**
     * The standard backbone letters among the search's declared series. Case is significant:
     * uppercase {@code Y} is MSFragger's peptide-remainder series, not the backbone {@code y}, and
     * is annotated as a peptide remainder rather than as a backbone ion. Custom names are resolved
     * separately, so anything else here is a series this result cannot generate.
     */
    private List<String> standardLetters(List<String> declared) {
        List<String> out = new ArrayList<>();
        List<String> unknown = new ArrayList<>();
        for (String s : declared) {
            if (s.length() == 1 && "abcxyz".contains(s)) {
                if (!out.contains(s)) out.add(s);
            } else if (s.equals("Y")) {
                continue; // handled as a peptide remainder, not a backbone series
            } else if (!hasCustom(s)) {
                unknown.add(s);
            }
        }
        if (!unknown.isEmpty()) {
            System.err.println("WARNING: fragment_ion_series names series this result does not "
                    + "declare and cannot generate: " + unknown);
        }
        return out;
    }

    private boolean hasCustom(String name) {
        for (CustomIon c : search.customIons) if (c.name.equals(name)) return true;
        return false;
    }

    /**
     * Load the glycan residue and modification databases, when FragPipe supplied them.
     * Pre-populates {@link FragmentAnnotator#glycoShortNames} for all loaded residues so label
     * generation is deterministic and thread-safe during annotation.
     */
    private void initGlycanResiduesMap() {
        if (glycanResiduesPath != null && !glycanResiduesPath.isEmpty()) {
            glycanResiduesMap = GlycanParser.parseGlycoResiduesDB(glycanResiduesPath);
        }
        if (glycanModsPath != null && !glycanModsPath.isEmpty()) {
            HashMap<String, umich.ms.glyco.GlycanMod> modsMap =
                    GlycanParser.parseGlycoModsDB(glycanModsPath, glycanResiduesMap.size(), glycanResiduesMap);
            glycanResiduesMap.putAll(modsMap);
        }
        if (glycanResiduesMap.isEmpty()) {
            if (search.nglycanMode) {
                System.err.println("WARNING: this is an N-glycan search but no glycan residue "
                        + "database was supplied, so glycan ions will be named by mass rather than "
                        + "by composition.");
            }
            return;
        }
        // Pre-register short names for all loaded residues before threads start.
        for (GlycanResidue residue : glycanResiduesMap.values()) {
            FragmentAnnotator.getOrCreateShortName(residue);
        }
        System.out.println("Loaded " + glycanResiduesMap.size() + " glycan residue/mod definitions.");
    }

    /** Per-category ion-match storage, in the order the columns are written. */
    private LinkedHashMap<IonCategory, ArrayList<IonMatch>[]> buildIonMatchesMap(int size) {
        LinkedHashMap<IonCategory, ArrayList<IonMatch>[]> map = new LinkedHashMap<>();
        for (IonCategory c : categories) map.put(c, new ArrayList[size]);
        return map;
    }

    private ScanCollectionDefault openMzmlFile(String spectrumName) {
        if (!resultProcessor.spectrumFileMap.containsKey(spectrumName)) return null;
        File eachFile = new File(resultProcessor.spectrumFileMap.get(spectrumName));
        if (eachFile.exists() && eachFile.getName().endsWith(".mzML")) {
            System.out.println("Reading mzML: " + spectrumName);
            MZMLFile mzmlFile = new MZMLFile(resultProcessor.spectrumFileMap.get(spectrumName));
            mzmlFile.setNumThreadsForParsing(threadNum);
            ScanCollectionDefault scans = new ScanCollectionDefault();
            scans.setDefaultStorageStrategy(StorageStrategy.SOFT);
            scans.isAutoloadSpectra(true);
            scans.setDataSource(mzmlFile);
            try {
                scans.loadData(LCMSDataSubset.STRUCTURE_ONLY);
                return scans;
            } catch (FileParsingException e) {
                e.printStackTrace();
            }
        }
        return null;
    }

    private ArrayList<String[]> readOnePSM(String expNum) throws IOException {
        System.out.println("Reading " + expNum);
        File onePSMTable = resultProcessor.resultsDict.get(expNum).get(1);
        ArrayList<String[]> onePSMData = new ArrayList<>();
        if (checkFileOpen(onePSMTable)) {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(onePSMTable));
            bufferedReader.readLine(); // skip header
            String line;
            while ((line = bufferedReader.readLine()) != null) {
                // Use limit -1 to preserve trailing empty fields, otherwise rows with
                // empty trailing columns (e.g. missing Protein Description / Mapped Genes)
                // would be written shorter than the header and shift the appended ion
                // columns to the left.
                onePSMData.add(line.split("\t", -1));
            }
            bufferedReader.close();
        }
        return onePSMData;
    }

    private HashMap<String, ArrayList<Integer>> groupPSMsByFile(ArrayList<String[]> onePSMData) {
        HashMap<String, ArrayList<Integer>> fileToIndices = new HashMap<>();
        for (int i = 0; i < onePSMData.size(); i++) {
            String fileName = onePSMData.get(i)[0].split("\\.")[0];
            fileToIndices.computeIfAbsent(fileName, k -> new ArrayList<>()).add(i);
        }
        return fileToIndices;
    }

    /** Load a spectrum's mz/intensity arrays from a scan collection. Returns null on failure. */
    private double[][] loadScanArrays(int scanNum, ScanCollectionDefault scans) throws FileParsingException {
        IScan iScan = scans.getScanByNum(scanNum);
        if (iScan == null) return null;
        ISpectrum iSpectrum = iScan.fetchSpectrum();
        if (iSpectrum == null) return null;
        double[] mzs = iSpectrum.getMZs();
        double[] ins  = iSpectrum.getIntensities();
        // specMzs must be sorted for binary search in FragmentAnnotator
        // batmass-io returns them sorted, but ensure:
        if (!isSorted(mzs)) sortParallel(mzs, ins);
        return new double[][]{mzs, ins};
    }

    private static boolean isSorted(double[] arr) {
        for (int i = 1; i < arr.length; i++) if (arr[i] < arr[i-1]) return false;
        return true;
    }

    private static void sortParallel(double[] mzs, double[] ins) {
        int n = mzs.length;
        Integer[] idx = new Integer[n];
        for (int i = 0; i < n; i++) idx[i] = i;
        Arrays.sort(idx, Comparator.comparingDouble(i -> mzs[i]));
        double[] tmpMz = mzs.clone(), tmpIn = ins.clone();
        for (int i = 0; i < n; i++) { mzs[i] = tmpMz[idx[i]]; ins[i] = tmpIn[idx[i]]; }
    }

    // Format one set of 5 ion annotation columns (tab-prefixed) for a single PSM row.
    private String formatIonColumns(ArrayList<IonMatch> matches, DecimalFormat df,
                                    DecimalFormat dfInt, DecimalFormat dfPpm) {
        if (matches == null || matches.isEmpty()) return "\t\t\t\t\t";
        ArrayList<String> ionsNames    = new ArrayList<>();
        ArrayList<String> ionsMz       = new ArrayList<>();
        ArrayList<String> ionsInt      = new ArrayList<>();
        ArrayList<String> ionsTheoMz   = new ArrayList<>();
        ArrayList<String> ionsPpmError = new ArrayList<>();
        for (IonMatch ionMatch : matches) {
            double ppmError = (ionMatch.theoMz - ionMatch.peakMz) / ionMatch.theoMz * 1e6;
            ionsNames.add(ionMatch.getPeakAnnotation());
            ionsMz.add(df.format(ionMatch.peakMz));
            ionsInt.add(dfInt.format(ionMatch.peakIntensity));
            ionsTheoMz.add(df.format(ionMatch.theoMz));
            ionsPpmError.add(dfPpm.format(ppmError));
        }
        return "\t" + ionsNames + "\t" + ionsMz + "\t" + ionsInt + "\t" + ionsTheoMz + "\t" + ionsPpmError;
    }

    /** The five column names of one category, optionally as the paired-scan copy. */
    private void appendHeader(StringBuilder sb, IonCategory category, boolean paired) {
        String prefix = (paired ? "paired_" : "") + category.prefix;
        sb.append("\t").append(prefix).append("ions")
          .append("\t").append(prefix).append("ion_mz")
          .append("\t").append(prefix).append("ion_int")
          .append("\t").append(prefix).append("ion_theo_mz")
          .append("\t").append(prefix).append("ion_ppm_error");
    }

    private void writeOneExperiment(String expNum, ArrayList<String[]> onePSMData,
                                    LinkedHashMap<IonCategory, ArrayList<IonMatch>[]> ionMatchesMap,
                                    LinkedHashMap<IonCategory, ArrayList<IonMatch>[]> pairedIonMatchesMap) {
        DecimalFormat df    = new DecimalFormat("#.####");
        DecimalFormat dfInt = new DecimalFormat("#.#");
        DecimalFormat dfPpm = new DecimalFormat("#.##");
        System.out.println("Writing " + expNum);
        File onePSMTable = resultProcessor.resultsDict.get(expNum).get(1);
        File onePSMTableWithMatch = new File(onePSMTable.getAbsolutePath().replace("psm.tsv", "psm_with_match.tsv"));
        try {
            BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(onePSMTableWithMatch));
            BufferedReader bufferedReader = new BufferedReader(new FileReader(onePSMTable));
            String line = bufferedReader.readLine();
            int columnNum = line.split("\t", -1).length;

            StringBuilder headerSuffix = new StringBuilder();
            for (IonCategory category : ionMatchesMap.keySet()) appendHeader(headerSuffix, category, false);
            if (pairedIonMatchesMap != null) {
                for (IonCategory category : pairedIonMatchesMap.keySet()) appendHeader(headerSuffix, category, true);
            }
            bufferedWriter.write(line.stripTrailing() + headerSuffix + "\n");
            bufferedReader.close();

            for (int i = 0; i < onePSMData.size(); i++) {
                String[] lineSplit = onePSMData.get(i);
                bufferedWriter.write(String.join("\t", lineSplit));
                for (int pad = lineSplit.length; pad < columnNum; pad++) bufferedWriter.write("\t");
                for (Map.Entry<IonCategory, ArrayList<IonMatch>[]> entry : ionMatchesMap.entrySet()) {
                    bufferedWriter.write(formatIonColumns(entry.getValue()[i], df, dfInt, dfPpm));
                }
                if (pairedIonMatchesMap != null) {
                    for (Map.Entry<IonCategory, ArrayList<IonMatch>[]> entry : pairedIonMatchesMap.entrySet()) {
                        bufferedWriter.write(formatIonColumns(entry.getValue()[i], df, dfInt, dfPpm));
                    }
                }
                bufferedWriter.write("\n");
            }
            bufferedWriter.close();
            onePSMTable.delete();
            onePSMTableWithMatch.renameTo(onePSMTable);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private Runnable getOneAnnotation(int psmIndexCount, ArrayList<String[]> onePSMData,
                                      LinkedHashMap<IonCategory, ArrayList<IonMatch>[]> ionMatchesMap,
                                      LinkedHashMap<IonCategory, ArrayList<IonMatch>[]> pairedIonMatchesMap,
                                      ScanCollectionDefault scans) {
        return () -> {
            try {
                String[] onePSM = onePSMData.get(psmIndexCount);
                int scanNum = Integer.parseInt(onePSM[0].split("\\.")[1]);
                int chargeValue = Integer.parseInt(onePSM[resultProcessor.chargeIndex]);

                double[][] primaryArrays = loadScanArrays(scanNum, scans);
                if (primaryArrays == null) return;
                double[] primaryMzs = primaryArrays[0];
                double[] primaryIns = primaryArrays[1];

                // Load paired spectrum if available
                double[] pairedMzs = null, pairedIns = null;
                if (pairedIonMatchesMap != null && resultProcessor.pairedScanNumIndex >= 0
                        && resultProcessor.pairedScanNumIndex < onePSM.length) {
                    String pairedScanStr = onePSM[resultProcessor.pairedScanNumIndex].trim();
                    if (!pairedScanStr.isEmpty()) {
                        int pairedScanNum = Integer.parseInt(pairedScanStr);
                        double[][] pairedArrays = loadScanArrays(pairedScanNum, scans);
                        if (pairedArrays != null) {
                            pairedMzs = pairedArrays[0];
                            pairedIns = pairedArrays[1];
                        }
                    }
                }

                String assignedMod     = onePSM[resultProcessor.assignenModIndex];
                String peptideSequence = onePSM[resultProcessor.peptideSequenceIndex];
                ArrayList<ModificationMatch> mods = parseModifications(assignedMod, peptideSequence);
                double deltaMass = readDeltaMass(onePSM);
                Glycan glycan = readGlycan(onePSM);

                EnumMap<IonCategory, ArrayList<IonMatch>> primary = FragmentAnnotator.annotate(
                        peptideSequence, mods, deltaMass, chargeValue,
                        primaryMzs, primaryIns,
                        backboneSeries, labileSeries, search, addNeutralLoss, glycan);
                store(ionMatchesMap, primary, psmIndexCount);

                if (pairedIonMatchesMap != null && pairedMzs != null) {
                    EnumMap<IonCategory, ArrayList<IonMatch>> paired = FragmentAnnotator.annotate(
                            peptideSequence, mods, deltaMass, chargeValue,
                            pairedMzs, pairedIns,
                            backboneSeries, labileSeries, search, addNeutralLoss, glycan);
                    store(pairedIonMatchesMap, paired, psmIndexCount);
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        };
    }

    private void store(LinkedHashMap<IonCategory, ArrayList<IonMatch>[]> target,
                       EnumMap<IonCategory, ArrayList<IonMatch>> matches, int psmIndex) {
        for (Map.Entry<IonCategory, ArrayList<IonMatch>[]> entry : target.entrySet()) {
            entry.getValue()[psmIndex] = matches.get(entry.getKey());
        }
    }

    /**
     * The PSM's unlocalized delta mass. FragPipe writes a localized mass offset into
     * Assigned Modifications and clears this, so a value here is a modification the search could
     * not place — it still produces diagnostic and peptide remainder ions, which need no position.
     */
    private double readDeltaMass(String[] onePSM) {
        int idx = resultProcessor.deltaMassIndex;
        if (idx < 0 || idx >= onePSM.length) return 0.0;
        Double d = SearchParams.parseDouble(onePSM[idx]);
        return d == null ? 0.0 : d;
    }

    /**
     * This PSM's glycan composition, or null. Its presence is what selects composition labels over
     * mass labels — see docs/adr/0001. A glyco search whose composition column is missing entirely
     * is warned about once, because that is the case where composition labels were expected.
     */
    private Glycan readGlycan(String[] onePSM) {
        if (glycanResiduesMap.isEmpty()) return null;
        int idx = resultProcessor.glycanCompositionIndex;
        if (idx < 0) {
            if (search.nglycanMode && warnedNoGlycanColumn.compareAndSet(false, true)) {
                System.err.println("WARNING: this is an N-glycan search but psm.tsv has no "
                        + "'Total Glycan Composition' column, so glycan ions will be named by mass. "
                        + "Run PTM-Shepherd's glycan assignment to get composition names.");
            }
            return null;
        }
        if (idx >= onePSM.length) return null;
        String glycanStr = onePSM[idx].trim();
        if (glycanStr.isEmpty()) return null;
        return GlycanParser.parseGlycanString(glycanStr, glycanResiduesMap);
    }

    /**
     * Parse FragPipe "Assigned Modifications" string into ModificationMatch list.
     * Format examples: "15M(15.9949)", "1S(79.9663)", "n(42.0106)", "N-term(229.1629)"
     */
    private ArrayList<ModificationMatch> parseModifications(String assignedMod, String peptideSequence) {
        ArrayList<ModificationMatch> mods = new ArrayList<>();
        if (assignedMod == null || assignedMod.trim().isEmpty()) return mods;

        for (String eachMod : assignedMod.split(",")) {
            eachMod = eachMod.trim();
            if (eachMod.contains(":") || !eachMod.contains("(")) continue;

            double modMass = Double.parseDouble(
                    eachMod.substring(eachMod.lastIndexOf('(') + 1, eachMod.lastIndexOf(')')));

            String ptmName;
            int position;

            if (eachMod.contains("n") && !eachMod.toLowerCase().contains("n-term")
                    && Character.isLowerCase(eachMod.charAt(0))) {
                // "n(42.0106)"
                ptmName  = modMass + " of N-term";
                position = 1;
            } else if (eachMod.toLowerCase().contains("n-term")) {
                // "N-term(229.1629)"
                String namedPtm = checkReporter(modMass, "N-term");
                ptmName  = namedPtm != null ? namedPtm : modMass + " of N-term";
                position = 1;
            } else if (eachMod.toLowerCase().contains("c-term")) {
                ptmName  = modMass + " of C-term";
                position = peptideSequence.length();
            } else {
                // e.g. "15M(15.9949)" → AA at position 15, residue 'M'
                String before = eachMod.substring(0, eachMod.lastIndexOf('('));
                String modAA  = before.substring(before.length() - 1);
                String posStr = before.substring(0, before.length() - 1).trim();
                position = Integer.parseInt(posStr);
                String namedPtm = checkReporter(modMass, modAA);
                ptmName  = namedPtm != null ? namedPtm : modMass + " of " + modAA;
            }

            mods.add(new ModificationMatch(ptmName, position, modMass));
        }
        return mods;
    }

    /** Map common reporter-tag masses to their conventional names. */
    private String checkReporter(double modMass, String modAA) {
        if (modAA.equals("N-term")) {
            if (Math.abs(modMass - 229.1629) <= 0.1) return "TMT 10-plex of peptide N-term";
            if (Math.abs(modMass - 144.1)    <= 0.1) return "iTRAQ 4-plex of peptide N-term";
        }
        return null;
    }

    private Boolean checkFileOpen(File eachFile) {
        return Files.exists(eachFile.toPath()) && Files.isRegularFile(eachFile.toPath())
                && Files.isReadable(eachFile.toPath());
    }
}
