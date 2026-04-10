import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scan.StorageStrategy;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.mzml.MZMLFile;

import java.io.*;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class ExportFragments {

    private int threadNum = 1;
    private ResultProcessor resultProcessor;
    private List<String> ionsTypeArray;
    private boolean hasNgly;
    private boolean hasOgly;

    public ExportFragments(File resultsFolder, int threadNum, List<String> ionsTypeArray) throws IOException {
        this.threadNum = threadNum;
        this.ionsTypeArray = ionsTypeArray;
        this.hasNgly = ionsTypeArray.contains("ngly") || ionsTypeArray.contains("gly");
        this.hasOgly = ionsTypeArray.contains("ogly") || ionsTypeArray.contains("gly");

        importData(resultsFolder);
    }

    private void importData(File resultsFolder) throws IOException {

        resultProcessor = new ResultProcessor(resultsFolder);

        if (resultProcessor.manifestFile != null) {

            if (resultProcessor.psmIndexToName.containsValue("ionint") && resultProcessor.psmIndexToName.containsValue("ionmz")){
                System.out.println("The file has already been annotated.");
                System.exit(1);
            }

            try {
                ExecutorService executorService = Executors.newFixedThreadPool(threadNum);

                for (String expNum : resultProcessor.resultsDict.keySet()) {
                    ArrayList<String[]> onePSMData = readOnePSM(expNum);
                    if (onePSMData == null || onePSMData.isEmpty()) continue;

                    HashMap<String, ArrayList<Integer>> fileToIndices = groupPSMsByFile(onePSMData);

                    LinkedHashMap<String, ArrayList<IonMatch>[]> ionMatchesMap = buildIonMatchesMap(onePSMData.size());
                    LinkedHashMap<String, ArrayList<IonMatch>[]> pairedIonMatchesMap =
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

    // Build per-type ionMatches storage.
    // ngly_ori – original glycan mass at N sites
    // ngly_0   – 0.0 of N (bare peptide backbone)
    // ngly_203 – 203.079 of N (core GlcNAc)
    // ogly_ori – original glycan mass at S/T sites
    // ogly_0   – 0.0 of ST (bare peptide backbone)
    private LinkedHashMap<String, ArrayList<IonMatch>[]> buildIonMatchesMap(int size) {
        LinkedHashMap<String, ArrayList<IonMatch>[]> map = new LinkedHashMap<>();
        if (hasNgly) {
            map.put("ngly_ori", new ArrayList[size]);
            map.put("ngly_0",   new ArrayList[size]);
            map.put("ngly_203", new ArrayList[size]);
        }
        if (hasOgly) {
            map.put("ogly_ori", new ArrayList[size]);
            map.put("ogly_0",   new ArrayList[size]);
        }
        if (!hasNgly && !hasOgly) map.put("default", new ArrayList[size]);
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

    private void writeOneExperiment(String expNum, ArrayList<String[]> onePSMData,
                                    LinkedHashMap<String, ArrayList<IonMatch>[]> ionMatchesMap,
                                    LinkedHashMap<String, ArrayList<IonMatch>[]> pairedIonMatchesMap) {
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
            for (String type : ionMatchesMap.keySet()) {
                String prefix = type.equals("default") ? "" : type + "_";
                headerSuffix.append("\t").append(prefix).append("ions")
                            .append("\t").append(prefix).append("ion_mz")
                            .append("\t").append(prefix).append("ion_int")
                            .append("\t").append(prefix).append("ion_theo_mz")
                            .append("\t").append(prefix).append("ion_ppm_error");
            }
            if (pairedIonMatchesMap != null) {
                for (String type : pairedIonMatchesMap.keySet()) {
                    String prefix = "paired_" + (type.equals("default") ? "" : type + "_");
                    headerSuffix.append("\t").append(prefix).append("ions")
                                .append("\t").append(prefix).append("ion_mz")
                                .append("\t").append(prefix).append("ion_int")
                                .append("\t").append(prefix).append("ion_theo_mz")
                                .append("\t").append(prefix).append("ion_ppm_error");
                }
            }
            bufferedWriter.write(line.stripTrailing() + headerSuffix + "\n");
            bufferedReader.close();

            for (int i = 0; i < onePSMData.size(); i++) {
                String[] lineSplit = onePSMData.get(i);
                bufferedWriter.write(String.join("\t", lineSplit));
                for (int pad = lineSplit.length; pad < columnNum; pad++) bufferedWriter.write("\t");
                for (Map.Entry<String, ArrayList<IonMatch>[]> entry : ionMatchesMap.entrySet()) {
                    bufferedWriter.write(formatIonColumns(entry.getValue()[i], df, dfInt, dfPpm));
                }
                if (pairedIonMatchesMap != null) {
                    for (Map.Entry<String, ArrayList<IonMatch>[]> entry : pairedIonMatchesMap.entrySet()) {
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
                                      LinkedHashMap<String, ArrayList<IonMatch>[]> ionMatchesMap,
                                      LinkedHashMap<String, ArrayList<IonMatch>[]> pairedIonMatchesMap,
                                      ScanCollectionDefault scans) {
        boolean anyGlycan = hasNgly || hasOgly;
        boolean neutral   = ionsTypeArray.contains("neu");
        return () -> {
            try {
                {
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

                    String assignedMod    = onePSM[resultProcessor.assignenModIndex];
                    String peptideSequence = onePSM[resultProcessor.peptideSequenceIndex];
                    ArrayList<ModificationMatch> originalMods = parseModifications(assignedMod, peptideSequence);

                    for (Map.Entry<String, ArrayList<IonMatch>[]> entry : ionMatchesMap.entrySet()) {
                        String type = entry.getKey();
                        ArrayList<IonMatch>[] typeIonMatches = entry.getValue();

                        ArrayList<ModificationMatch> activeMods = buildModsForType(originalMods, peptideSequence, type);

                        boolean addGlycan = anyGlycan && (type.equals("ngly_ori") || type.equals("ogly_ori"));

                        typeIonMatches[psmIndexCount] = FragmentAnnotator.annotate(
                                peptideSequence, activeMods, chargeValue,
                                primaryMzs, primaryIns,
                                ionsTypeArray, addGlycan, neutral);

                        if (pairedIonMatchesMap != null) {
                            ArrayList<IonMatch>[] pairedTypeIonMatches = pairedIonMatchesMap.get(type);
                            if (pairedMzs != null) {
                                pairedTypeIonMatches[psmIndexCount] = FragmentAnnotator.annotate(
                                        peptideSequence, activeMods, chargeValue,
                                        pairedMzs, pairedIns,
                                        ionsTypeArray, addGlycan, neutral);
                            }
                        }
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        };
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

    /**
     * Build the mod list for a specific annotation type variant.
     * For N-glycan types, replaces glycan mods (N with mass>500) with the
     * appropriate variant mass. For O-glycan types, similar for S/T>100.
     */
    private ArrayList<ModificationMatch> buildModsForType(
            ArrayList<ModificationMatch> originalMods, String peptideSequence, String type) {

        if (type.equals("default")) return originalMods;

        ArrayList<ModificationMatch> newMods = new ArrayList<>();
        for (ModificationMatch mm : originalMods) {
            String ptm = mm.getTheoreticPtm();
            if (!ptm.contains(" of ")) { newMods.add(mm); continue; }
            String residue = ptm.split(" of ")[1];
            double mass    = mm.getMass();
            int site       = mm.getModificationSite();

            if (type.startsWith("ngly") && residue.equals("N")
                    && (mass > 500 || ptm.startsWith("203.079") || ptm.startsWith("0.0"))) {
                switch (type) {
                    case "ngly_ori": newMods.add(mm); break;
                    case "ngly_0":
                        newMods.add(new ModificationMatch("0.0 of N", site, 0.0)); break;
                    case "ngly_203":
                        newMods.add(new ModificationMatch("203.079 of N", site, 203.07937)); break;
                }
            } else if (type.startsWith("ogly")
                    && (residue.equals("S") || residue.equals("T") || residue.equals("ST"))
                    && (mass > 100 || ptm.startsWith("0.0"))) {
                switch (type) {
                    case "ogly_ori": newMods.add(mm); break;
                    case "ogly_0":
                        newMods.add(new ModificationMatch("0.0 of ST", site, 0.0)); break;
                }
            } else {
                newMods.add(mm);
            }
        }
        return newMods;
    }


private Boolean checkFileOpen(File eachFile) {
        return Files.exists(eachFile.toPath()) && Files.isRegularFile(eachFile.toPath())
                && Files.isReadable(eachFile.toPath());
    }
}
