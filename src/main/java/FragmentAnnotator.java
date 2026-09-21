import umich.ms.glyco.Glycan;
import umich.ms.glyco.GlycanResidue;

import java.util.*;

/**
 * Fragment ion annotator: calculates theoretical ions and matches them to observed spectrum peaks
 * within the search's own fragment tolerance.
 *
 * <p>All masses are monoisotopic (Da).
 * <pre>
 *   b neutral = sum(N-term residues + mods)
 *   y neutral = sum(C-term residues + mods) + H2O
 *   a = b - CO,  c = b + NH3,  x = y + CO2 - H2O,  z-radical = y - NH2
 * </pre>
 *
 * <p>Every series — standard or user-defined — is a {@link Series}: a label, a terminus and a
 * shift from that terminus' ordinary ion. No code may test an ion type against a literal list of
 * letters, because a custom series named {@code zOne} passes no such test.
 */
public class FragmentAnnotator {

    // Physical constants ----------------------------------------------------
    public static final double PROTON    = 1.007276466812; // proton (not H atom)
    public static final double H_ATOM    = 1.00782503207;  // H atom (proton + electron)
    public static final double H2O       = 18.01056468326; // 2H + O
    public static final double NH3       = 17.02654910101; // N + 3H
    public static final double H3PO4     = 97.97689540;
    public static final double CO        = 27.99491461956; // C + O
    public static final double CO2       = 43.98982923912; // C + 2O
    // z-radical offset from y: y - NH2 = y - NH3 + H_atom
    private static final double Z_DOT_DELTA = NH3 - H_ATOM; // = N + 2H = 16.01872406894

    // Amino acid residue masses ---------------------------------------------
    private static final Map<Character, Double> AA_MASS = new HashMap<>();
    static {
        AA_MASS.put('A', 71.03711);
        AA_MASS.put('R', 156.10111);
        AA_MASS.put('N', 114.04293);
        AA_MASS.put('D', 115.02694);
        AA_MASS.put('C', 103.00919);
        AA_MASS.put('E', 129.04259);
        AA_MASS.put('Q', 128.05858);
        AA_MASS.put('G', 57.02146);
        AA_MASS.put('H', 137.05891);
        AA_MASS.put('I', 113.08406);
        AA_MASS.put('L', 113.08406);
        AA_MASS.put('K', 128.09496);
        AA_MASS.put('M', 131.04049);
        AA_MASS.put('F', 147.06841);
        AA_MASS.put('P', 97.05276);
        AA_MASS.put('S', 87.03203);
        AA_MASS.put('T', 101.04768);
        AA_MASS.put('W', 186.07931);
        AA_MASS.put('Y', 163.06333);
        AA_MASS.put('V', 99.06841);
        AA_MASS.put('U', 150.95363); // selenocysteine
        AA_MASS.put('O', 237.14773); // pyrrolysine
    }

    // Glycan monosaccharide residue masses ----------------------------------
    // Residue mass = intact monosaccharide mass - H2O
    public static final double HEXNAC_RESIDUE = 203.07937; // HexNAc (GlcNAc/GalNAc)
    public static final double HEX_RESIDUE    = 162.05282; // Hex (Man/Gal/Glc)
    public static final double NEUAC_RESIDUE  = 291.09542; // NeuAc (Sialic acid)
    public static final double NEUGC_RESIDUE  = 307.09033; // NeuGc
    public static final double DHEX_RESIDUE   = 146.05791; // dHex (Fucose)

    /**
     * Named singly-charged oxonium ions: {label, m/z}, where m/z = sum(residue masses) + PROTON.
     *
     * <p>On a glycopeptide PSM these are merged with the m/z values the search declared. The
     * search's own list is composition-specific and so omits legitimate markers — an offset for a
     * HexNAc-only glycan declares no sialic-acid oxonium — while this table names ions the
     * parameter leaves anonymous. Neither alone is enough, so the union is annotated and this
     * table supplies the names.
     */
    private static final Object[][] OXONIUM_IONS = {
        {"B_N1",         HEXNAC_RESIDUE + PROTON},
        {"B_N1-H2O",     HEXNAC_RESIDUE + PROTON - H2O},
        {"B_N1-2H2O",    HEXNAC_RESIDUE + PROTON - 2 * H2O},
        {"B_H1",         HEX_RESIDUE + PROTON},
        {"B_A1",         NEUAC_RESIDUE + PROTON},
        {"B_A1-H2O",     NEUAC_RESIDUE + PROTON - H2O},
        {"B_A1-2H2O",    NEUAC_RESIDUE + PROTON - 2 * H2O},
        {"B_G1",         NEUGC_RESIDUE + PROTON},
        {"B_F1",         DHEX_RESIDUE + PROTON},
        {"B_N2",         2 * HEXNAC_RESIDUE + PROTON},
        {"B_N2-H2O",     2 * HEXNAC_RESIDUE + PROTON - H2O},
        {"B_N1H1",       HEX_RESIDUE + HEXNAC_RESIDUE + PROTON},
        {"B_N1H1-H2O",   HEX_RESIDUE + HEXNAC_RESIDUE + PROTON - H2O},
        {"B_N1A1",       HEXNAC_RESIDUE + NEUAC_RESIDUE + PROTON},
        {"B_N2H1",       2 * HEXNAC_RESIDUE + HEX_RESIDUE + PROTON},
        {"B_N2H2",       2 * HEXNAC_RESIDUE + 2 * HEX_RESIDUE + PROTON},
        {"B_N2H3",       2 * HEXNAC_RESIDUE + 3 * HEX_RESIDUE + PROTON},
        {"B_N1F1",       DHEX_RESIDUE + HEXNAC_RESIDUE + PROTON},
    };

    // Glycan residue short names for compact Y ion labels.
    public static final HashMap<String, String> glycoShortNames = new HashMap<>();
    static {
        glycoShortNames.put("HexNAc", "N");
        glycoShortNames.put("Hex",    "H");
        glycoShortNames.put("Fuc",    "F");
        glycoShortNames.put("NeuAc",  "A");
        glycoShortNames.put("NeuGc",  "G");
    }

    // Peak matching tolerance -----------------------------------------------
    // Set once from the search parameters before any worker thread starts, and read-only after.
    private static volatile double tolValue = SearchParams.DEFAULT_TOL_PPM;
    private static volatile boolean tolIsDa = false;

    /**
     * Sets the peak-matching tolerance from the search's own parameters. Must be called before
     * worker threads are spawned; the fields are read-only afterwards.
     */
    public static void configureTolerance(double value, boolean isDa) {
        tolValue = value;
        tolIsDa = isDa;
    }

    // Series ----------------------------------------------------------------

    /**
     * One ion series: a label, a terminus, and the shift from that terminus' ordinary ion (b for
     * N-terminal, y for C-terminal). Standard and user-defined series differ only in their values.
     */
    public static class Series {
        public final String label;
        public final boolean nterm;
        public final double shift;

        public Series(String label, boolean nterm, double shift) {
            this.label = label;
            this.nterm = nterm;
            this.shift = shift;
        }
    }

    /** The standard series by their MSFragger letters. */
    private static Series standard(String letter) {
        switch (letter) {
            case "a": return new Series("a", true,  -CO);
            case "b": return new Series("b", true,  0.0);
            case "c": return new Series("c", true,  NH3);
            case "x": return new Series("x", false, CO2 - H2O);
            case "y": return new Series("y", false, 0.0);
            case "z": return new Series("z", false, -Z_DOT_DELTA);
            default:  return null;
        }
    }

    /**
     * The generator series for a set of labels: a standard letter, or a custom series by name.
     * An unknown label is dropped with a warning — a series this result never declared has no
     * offset to generate ions from, and guessing one would annotate peaks under a false name.
     */
    public static List<Series> resolveSeries(List<String> labels, List<CustomIon> customIons) {
        List<Series> out = new ArrayList<>();
        for (String label : labels) {
            Series s = standard(label);
            if (s != null) {
                out.add(s);
                continue;
            }
            CustomIon ci = findCustom(customIons, label);
            if (ci != null) out.add(new Series(ci.name, ci.nterm, ci.shiftFromBase()));
        }
        return out;
    }

    private static CustomIon findCustom(List<CustomIon> customIons, String label) {
        if (customIons == null) return null;
        for (CustomIon c : customIons) if (c.name.equals(label)) return c;
        return null;
    }

    // Annotation ------------------------------------------------------------

    /**
     * A labile modification found on this PSM, paired with the offset that explains it.
     *
     * <p>{@code site} is the 1-based residue position, or 0 for an unlocalized modification, which
     * has no position and so can only produce the two non-positional ion kinds. {@code modMass} is
     * what must be removed from the peptide to strip it — zero when unlocalized, because such a
     * mass never entered the computed peptide mass in the first place.
     */
    private static class Hit {
        final int site;
        final double modMass;
        final LabileOffset off;

        Hit(int site, double modMass, LabileOffset off) {
            this.site = site;
            this.modMass = modMass;
            this.off = off;
        }
    }

    /**
     * Annotate a spectrum, sorting the matches into the four ion categories.
     *
     * @param sequence         peptide sequence (single-letter codes, upper-case)
     * @param mods             modifications (position 1-indexed)
     * @param deltaMass        the PSM's unlocalized delta mass, or 0
     * @param precursorCharge  precursor charge state
     * @param specMzs          sorted observed m/z values
     * @param specInts         observed intensities (same order as specMzs)
     * @param backbone         backbone series to generate (standard and custom)
     * @param labileBases      series that carry fragment remainders
     * @param search           the search's parameters; labile ions are skipped when not labile
     * @param addNeutralLoss   whether to annotate neutral-loss variants of b and y
     * @param glycan           this PSM's glycan composition, or null; its presence is what selects
     *                         composition labels over mass labels
     */
    public static EnumMap<IonCategory, ArrayList<IonMatch>> annotate(
            String sequence,
            ArrayList<ModificationMatch> mods,
            double deltaMass,
            int precursorCharge,
            double[] specMzs,
            double[] specInts,
            List<Series> backbone,
            List<Series> labileBases,
            SearchParams search,
            boolean addNeutralLoss,
            Glycan glycan) {

        EnumMap<IonCategory, ArrayList<IonMatch>> result = new EnumMap<>(IonCategory.class);
        for (IonCategory c : IonCategory.values()) result.put(c, new ArrayList<>());
        if (specMzs == null || specMzs.length == 0) return result;

        int n = sequence.length();

        // Per-position mod mass arrays (1-indexed: index 1..n for residues)
        double nTermDelta = 0.0;
        double cTermDelta = 0.0;
        double[] resDelta = new double[n + 1];

        for (ModificationMatch mm : mods) {
            double mass = mm.getMass();
            if (mm.isNTerminal()) {
                nTermDelta += mass;
            } else if (mm.isCTerminal()) {
                cTermDelta += mass;
            } else {
                int pos = mm.getModificationSite();
                if (pos >= 1 && pos <= n) resDelta[pos] += mass;
            }
        }

        // Prefix sums including residue + mod mass at each position
        double[] prefix = new double[n + 1];
        for (int i = 1; i <= n; i++) {
            double aaMass = AA_MASS.getOrDefault(sequence.charAt(i - 1), 0.0);
            prefix[i] = prefix[i - 1] + aaMass + resDelta[i];
        }

        double totalPeptideMass = prefix[n] + H2O + nTermDelta + cTermDelta;
        int maxCharge = Math.min(Math.max(precursorCharge, 1), 4);

        // Backbone ions -----------------------------------------------------
        List<double[]> neutralLosses = new ArrayList<>();
        if (addNeutralLoss) {
            neutralLosses.add(new double[]{H2O, 0});
            neutralLosses.add(new double[]{NH3, 1});
            for (ModificationMatch mm : mods) {
                if (isPhospho(mm)) {
                    neutralLosses.add(new double[]{H3PO4, 2});
                    break;
                }
            }
        }
        String[] lossNames = {"-H2O", "-NH3", "-H3PO4"};

        ArrayList<IonMatch> backboneOut = result.get(IonCategory.BACKBONE);
        for (int i = 1; i < n; i++) {
            double bNeutral = prefix[i] + nTermDelta;
            double yNeutral = prefix[n] - prefix[i] + H2O + cTermDelta;

            for (Series s : backbone) {
                double neutral = (s.nterm ? bNeutral : yNeutral) + s.shift;
                int index = s.nterm ? i : n - i;
                for (int z = 1; z <= maxCharge; z++) {
                    String name = s.label + index + chargeStr(z);
                    addMatch(backboneOut, specMzs, specInts, name, (neutral + z * PROTON) / z);

                    // Neutral losses stay on b and y only, as before.
                    if (!s.label.equals("b") && !s.label.equals("y")) continue;
                    for (double[] loss : neutralLosses) {
                        addMatch(backboneOut, specMzs, specInts,
                                name + lossNames[(int) loss[1]],
                                (neutral - loss[0] + z * PROTON) / z);
                    }
                }
            }
        }

        // A search with no per-mass offsets may still have global ion lists a glycan can use: a
        // glyco search keeps its glycan masses in mass_offsets_detailed and writes mass_offsets=0,
        // so with detailed offsets off there is nothing to match by mass and everything to annotate.
        if (!search.labile || (search.offsets.isEmpty() && search.globalOffset == null)) return result;

        // Which declared offsets this PSM actually carries ------------------
        List<Hit> hits = new ArrayList<>();
        for (ModificationMatch mm : mods) {
            int pos = mm.getModificationSite();
            char aa = (pos >= 1 && pos <= n) ? sequence.charAt(pos - 1) : 0;
            LabileOffset off = search.offsetFor(mm.getMass(), aa);
            if (off != null) hits.add(new Hit(pos, mm.getMass(), off));
        }
        // An offset the search could not localize survives as a bare delta mass. It has no
        // position, so it produces the two non-positional kinds only, and nothing to strip from
        // the peptide, because its mass never entered the computed peptide mass.
        if (Math.abs(deltaMass) > SearchParams.OFFSET_TOL) {
            LabileOffset off = search.offsetFor(deltaMass, (char) 0);
            if (off != null) hits.add(new Hit(0, 0.0, off));
        }
        // A glycan whose mass has no entry of its own still fragments. A glyco search keeps its
        // glycan masses in mass_offsets_detailed and leaves mass_offsets=0, so with detailed
        // offsets switched off there is no per-mass entry to find; the search's global ion lists
        // are what it used, and without this the glycopeptide would get no labile ions at all.
        if (glycan != null && search.globalOffset != null && !covers(hits, glycan.mass)) {
            for (ModificationMatch mm : mods) {
                if (Math.abs(mm.getMass() - glycan.mass) < SearchParams.OFFSET_TOL) {
                    hits.add(new Hit(mm.getModificationSite(), mm.getMass(), search.globalOffset));
                    break;
                }
            }
        }
        if (hits.isEmpty()) return result;

        annotateDiagnostic(result.get(IonCategory.DIAGNOSTIC), specMzs, specInts, hits, glycan != null);
        annotatePepRemainders(result.get(IonCategory.PEP_REMAINDER), specMzs, specInts,
                hits, totalPeptideMass, maxCharge, glycan);
        annotateFragRemainders(result.get(IonCategory.FRAG_REMAINDER), specMzs, specInts,
                hits, prefix, nTermDelta, cTermDelta, n, labileBases, maxCharge);
        return result;
    }

    /** Whether some hit already explains a modification of this mass. */
    private static boolean covers(List<Hit> hits, double mass) {
        for (Hit h : hits) if (Math.abs(h.modMass - mass) < SearchParams.OFFSET_TOL) return true;
        return false;
    }

    /**
     * Diagnostic ions, at the m/z the parameter states: it carries m/z, not neutral masses, so
     * these are used as given rather than protonated. Deduplicated across the modifications that
     * produced them — two residues carrying the same offset describe the same ion, and annotating
     * it twice would put two labels on one peak.
     */
    private static void annotateDiagnostic(ArrayList<IonMatch> out, double[] specMzs, double[] specInts,
                                           List<Hit> hits, boolean glycopeptide) {
        List<Double> mzs = new ArrayList<>();
        List<String> names = new ArrayList<>();
        for (Hit h : hits) {
            for (double mz : h.off.diagnostic) {
                // Snap a declared m/z onto the table entry that names it, BEFORE deduplicating.
                // Naming and deduplication must use one window: a wider naming window would give
                // the same B_ name to two entries the narrower dedup kept apart, and at these
                // tolerances both could then match one peak and label it twice.
                Object[] named = glycopeptide ? oxoniumEntry(mz) : null;
                if (named != null) addUniqueMz(mzs, names, (double) named[1], (String) named[0]);
                else addUniqueMz(mzs, names, mz, null);
            }
        }
        if (glycopeptide) {
            for (Object[] entry : OXONIUM_IONS) addUniqueMz(mzs, names, (double) entry[1], (String) entry[0]);
        }
        for (int i = 0; i < mzs.size(); i++) {
            String name = names.get(i);
            if (name == null) name = String.format(Locale.ROOT, "%.2f", mzs.get(i));
            addMatch(out, specMzs, specInts, name, mzs.get(i));
        }
    }

    /** Adds an m/z unless an equal one is already present, keeping the first name offered for it. */
    private static void addUniqueMz(List<Double> mzs, List<String> names, double mz, String name) {
        for (int i = 0; i < mzs.size(); i++) {
            if (Math.abs(mzs.get(i) - mz) < 1e-4) {
                if (names.get(i) == null && name != null) names.set(i, name);
                return;
            }
        }
        mzs.add(mz);
        names.add(name);
    }

    /**
     * The oxonium table entry naming this m/z, or null. A declared ion the table does not know
     * keeps its own m/z as its label, because the parameter names none of its ions and inventing
     * a name would claim a chemistry nothing recorded.
     */
    private static Object[] oxoniumEntry(double mz) {
        for (Object[] entry : OXONIUM_IONS) {
            if (Math.abs((double) entry[1] - mz) < OXONIUM_NAME_TOL) return entry;
        }
        return null;
    }

    /**
     * How far a declared diagnostic m/z may sit from a named oxonium ion and still be it. Wide
     * enough for the rounding MSFragger writes (its values agree with the table to about 1e-5),
     * narrow enough not to claim a neighbouring ion's name.
     */
    private static final double OXONIUM_NAME_TOL = 1e-3;

    /**
     * Peptide remainder ions: the intact peptide with its labile modifications stripped and a
     * remainder put back. Not positional, so one set per PSM.
     *
     * <p>On a glycopeptide PSM the glycan's own sub-compositions replace the declared remainder
     * masses, because they carry the same masses under names an analyst can read.
     */
    private static void annotatePepRemainders(ArrayList<IonMatch> out, double[] specMzs, double[] specInts,
                                              List<Hit> hits, double totalPeptideMass,
                                              int maxCharge, Glycan glycan) {
        if (glycan != null) {
            generateGlycanYIons(out, specMzs, specInts, totalPeptideMass - glycan.mass, glycan, maxCharge);
            return;
        }
        double stripped = totalPeptideMass;
        List<double[]> lists = new ArrayList<>();
        for (Hit h : hits) {
            if (h.off.peptide.length == 0) continue;
            stripped -= h.modMass;
            lists.add(h.off.peptide);
        }
        if (lists.isEmpty()) return;
        for (double sum : uniqueSums(lists)) {
            String label = "pep" + remainderTag(sum);
            for (int z = 1; z <= maxCharge; z++) {
                addMatch(out, specMzs, specInts, label + chargeStr(z), (stripped + sum + z * PROTON) / z);
            }
        }
    }

    /**
     * Fragment remainder ions: backbone ions whose labile modifications were partly or wholly
     * lost. Only a fragment that CONTAINS a labile modification can carry its remainder — a
     * fragment with none is an ordinary backbone ion and must not be annotated twice.
     *
     * <p>Where a fragment contains several labile modifications the remainders combine, and the
     * ion is labelled by their sum: the ion's identity is its mass, so two combinations summing
     * alike are one peak and take one label.
     */
    private static void annotateFragRemainders(ArrayList<IonMatch> out, double[] specMzs, double[] specInts,
                                               List<Hit> hits, double[] prefix,
                                               double nTermDelta, double cTermDelta, int n,
                                               List<Series> bases, int maxCharge) {
        if (bases.isEmpty()) return;
        List<Hit> positional = new ArrayList<>();
        for (Hit h : hits) if (h.site >= 1 && h.off.fragment.length > 0) positional.add(h);
        if (positional.isEmpty()) return;

        // The mass adjustment is the same for every series and charge at a given cleavage, so it
        // is computed once per cleavage per terminus rather than once per ion.
        double[][] nDelta = new double[n + 1][];
        double[][] cDelta = new double[n + 1][];
        String[][] nTag = new String[n + 1][];
        String[][] cTag = new String[n + 1][];
        for (int i = 1; i < n; i++) {
            // An N-terminal fragment of length i contains sites 1..i; a C-terminal one of length
            // n-i contains sites i+1..n.
            fillRemainders(positional, i, true, nDelta, nTag, i);
            fillRemainders(positional, i, false, cDelta, cTag, i);
        }

        for (int i = 1; i < n; i++) {
            double bNeutral = prefix[i] + nTermDelta;
            double yNeutral = prefix[n] - prefix[i] + H2O + cTermDelta;
            for (Series s : bases) {
                double neutral = (s.nterm ? bNeutral : yNeutral) + s.shift;
                int index = s.nterm ? i : n - i;
                double[] deltas = s.nterm ? nDelta[i] : cDelta[i];
                String[] tags = s.nterm ? nTag[i] : cTag[i];
                if (deltas == null) continue;
                for (int k = 0; k < deltas.length; k++) {
                    for (int z = 1; z <= maxCharge; z++) {
                        addMatch(out, specMzs, specInts, s.label + index + tags[k] + chargeStr(z),
                                (neutral + deltas[k] + z * PROTON) / z);
                    }
                }
            }
        }
    }

    /** Computes the distinct mass adjustments and labels for one cleavage at one terminus. */
    private static void fillRemainders(List<Hit> positional, int cleavage, boolean nterm,
                                       double[][] deltas, String[][] tags, int slot) {
        double contained = 0.0;
        List<double[]> lists = new ArrayList<>();
        for (Hit h : positional) {
            boolean inside = nterm ? h.site <= cleavage : h.site > cleavage;
            if (!inside) continue;
            contained += h.modMass;
            lists.add(h.off.fragment);
        }
        if (lists.isEmpty()) return;
        double[] sums = uniqueSums(lists);
        double[] d = new double[sums.length];
        String[] t = new String[sums.length];
        int n = 0;
        for (double sum : sums) {
            double delta = sum - contained;
            // A remainder that restores the full modification mass IS the backbone ion.
            if (Math.abs(delta) < 1e-6) continue;
            d[n] = delta;
            t[n] = remainderTag(sum);
            n++;
        }
        deltas[slot] = Arrays.copyOf(d, n);
        tags[slot] = Arrays.copyOf(t, n);
    }

    /**
     * The distinct sums of one value drawn from each list — the cross-product of several labile
     * modifications' remainders, collapsed to the masses it can actually produce.
     */
    static double[] uniqueSums(List<double[]> lists) {
        double[] sums = {0.0};
        for (double[] list : lists) {
            double[] next = new double[sums.length * Math.max(list.length, 1)];
            int n = 0;
            for (double s : sums) {
                for (double v : list) {
                    double sum = s + v;
                    boolean seen = false;
                    for (int i = 0; i < n; i++) {
                        if (Math.abs(next[i] - sum) < 1e-6) { seen = true; break; }
                    }
                    if (!seen) next[n++] = sum;
                }
            }
            sums = Arrays.copyOf(next, n);
        }
        return sums;
    }

    /**
     * A remainder mass as it appears on a label: signed and rounded to the nearest integer
     * ({@code "+203"}, {@code "-42"}, {@code "+0"}). Rounded because the exact mass is already in
     * the theoretical-m/z column.
     */
    static String remainderTag(double remainder) {
        long r = Math.round(remainder);
        return (r < 0 ? "-" : "+") + Math.abs(r);
    }

    /**
     * Generate glycan Y ions by enumerating every valid fragment sub-composition of the supplied
     * glycan. Only compositions that are actual subsets of the intact glycan are produced.
     */
    private static void generateGlycanYIons(
            ArrayList<IonMatch> result,
            double[] specMzs, double[] specInts,
            double barePeptideMass,
            Glycan glycan,
            int maxCharge) {

        List<LinkedHashMap<GlycanResidue, Integer>> fragComps = glycan.generateFragmentCompositions();
        for (LinkedHashMap<GlycanResidue, Integer> fragComp : fragComps) {
            double comboMass = Glycan.computeCompositionMass(fragComp);
            double yMass = barePeptideMass + comboMass;
            // Build label lazily — only when a peak match is found.
            String label = null;
            for (int z = 1; z <= maxCharge; z++) {
                double mz = (yMass + z * PROTON) / z;
                int peakIdx = findBestPeak(specMzs, mz);
                if (peakIdx >= 0) {
                    if (label == null) label = buildGlycanYLabel(fragComp);
                    result.add(new IonMatch(label + chargeStr(z), specMzs[peakIdx], specInts[peakIdx], mz));
                }
            }
        }
    }

    /**
     * Build a compact Y_ label from a fragment composition map, e.g. {HexNAc-2, Hex-1} to
     * {@code "Y_N2H1"}. The empty composition is the bare peptide, and is named {@code "Y_0"}
     * rather than left as a bare prefix.
     */
    private static String buildGlycanYLabel(LinkedHashMap<GlycanResidue, Integer> fragComp) {
        StringBuilder sb = new StringBuilder("Y_");
        for (Map.Entry<GlycanResidue, Integer> entry : fragComp.entrySet()) {
            if (entry.getValue() == null || entry.getValue() == 0) continue;
            sb.append(getOrCreateShortName(entry.getKey())).append(entry.getValue());
        }
        if (sb.length() == 2) sb.append('0');
        return sb.toString();
    }

    /**
     * Return the short name for a glycan residue, registering a new one if needed. Finds the
     * shortest prefix of {@code residue.name} that is not already in use. Synchronized to be safe
     * when called from multiple threads (residues are pre-registered at startup, so contention is
     * rare).
     */
    static synchronized String getOrCreateShortName(GlycanResidue residue) {
        String name = residue.name;
        String existing = glycoShortNames.get(name);
        if (existing != null) return existing;
        Set<String> usedValues = new HashSet<>(glycoShortNames.values());
        String shortName = name; // fallback: use full name if all prefixes are taken
        for (int len = 1; len <= name.length(); len++) {
            String candidate = name.substring(0, len);
            if (!usedValues.contains(candidate)) {
                shortName = candidate;
                break;
            }
        }
        glycoShortNames.put(name, shortName);
        return shortName;
    }

    // Precomputed charge suffix strings to avoid repeated String.repeat() allocations
    private static final String[] CHARGE_STR = {"", "", "++", "+++", "++++"};

    /** Try to match a theoretical m/z to the nearest observed peak within the search tolerance. */
    private static void addMatch(
            ArrayList<IonMatch> result,
            double[] specMzs, double[] specInts,
            String label, double theoMz) {

        if (theoMz <= 0) return;
        int best = findBestPeak(specMzs, theoMz);
        if (best >= 0) result.add(new IonMatch(label, specMzs[best], specInts[best], theoMz));
    }

    /**
     * Returns the index of the nearest peak within the search's fragment tolerance, or -1 if none.
     * The tolerance is the one MSFragger used, in ppm or in Da: a Da tolerance is an absolute
     * window and cannot be expressed as a ppm one.
     */
    private static int findBestPeak(double[] specMzs, double theoMz) {
        if (theoMz <= 0 || specMzs.length == 0) return -1;
        int lo = 0, hi = specMzs.length - 1;
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (specMzs[mid] < theoMz) lo = mid + 1;
            else hi = mid;
        }
        int best = lo;
        if (lo > 0 && Math.abs(specMzs[lo - 1] - theoMz) < Math.abs(specMzs[lo] - theoMz)) {
            best = lo - 1;
        }
        double diff = Math.abs(specMzs[best] - theoMz);
        if (tolIsDa) return diff <= tolValue ? best : -1;
        return diff / theoMz * 1e6 <= tolValue ? best : -1;
    }

    /** Returns charge indicator string: "" for z=1, "++" for z=2, etc. */
    private static String chargeStr(int z) {
        return z < CHARGE_STR.length ? CHARGE_STR[z] : repeat(z);
    }

    private static String repeat(int z) {
        StringBuilder sb = new StringBuilder(z);
        for (int i = 0; i < z; i++) sb.append('+');
        return sb.toString();
    }

    /** True when the modification represents a phosphorylation on S or T (~80 Da). */
    private static boolean isPhospho(ModificationMatch mm) {
        String ptm = mm.getTheoreticPtm();
        if (!ptm.contains(" of ")) return false;
        String residue = ptm.split(" of ")[1];
        if (!residue.equals("S") && !residue.equals("T")) return false;
        double mass = mm.getMass();
        return mass > 79.9 && mass < 80.01;
    }
}
