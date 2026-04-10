import java.util.*;

/**
 * Fragment ion annotator: calculates theoretical a/b/c/x/y/z• ions and glycan B/Y ions,
 * matches them to observed spectrum peaks within a 20 ppm window.
 *
 * All masses are monoisotopic (Da).
 *   b neutral = sum(N-term residues + mods)
 *   y neutral = sum(C-term residues + mods) + H2O
 *   a = b - CO,  c = b + NH3,  x = y + CO2 - H2O,  z• = y - NH2• (radical)
 */
public class FragmentAnnotator {

    // ── Physical constants ──────────────────────────────────────────────────
    public static final double PROTON    = 1.007276466812; // proton (not H atom)
    public static final double H_ATOM    = 1.00782503207;  // H atom (proton + electron)
    public static final double H2O       = 18.01056468326; // 2H + O
    public static final double NH3       = 17.02654910101; // N + 3H
    public static final double H3PO4     = 97.97689540;
    public static final double CO        = 27.99491461956; // C + O
    public static final double CO2       = 43.98982923912; // C + 2O
    // z• (radical) offset from y: y - NH2• = y - NH3 + H_atom
    // Equivalent: rewindMass - N_ATOM where N_ATOM = 14.0030740048
    private static final double Z_DOT_DELTA = NH3 - H_ATOM; // = N + 2H = 16.01872406894

    // ── Amino acid residue masses ───────────────────────────────────────────
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

    // ── Glycan monosaccharide residue masses ────────────────────────────────
    // Residue mass = intact monosaccharide mass - H2O
    public static final double HEXNAC_RESIDUE = 203.07937; // HexNAc (GlcNAc/GalNAc)
    public static final double HEX_RESIDUE    = 162.05282; // Hex (Man/Gal/Glc)
    public static final double NEUAC_RESIDUE  = 291.09542; // NeuAc (Sialic acid)
    public static final double NEUGC_RESIDUE  = 307.09033; // NeuGc
    public static final double DHEX_RESIDUE   = 146.05791; // dHex (Fucose)

    // Sugar single-letter codes (N=HexNAc, H=Hex, A=NeuAc, G=NeuGc, F=Fucose)
    private static final String[] SUGAR_LETTERS = {"N", "H", "A", "G", "F"};
    private static final double[] SUGAR_MASSES  = {HEXNAC_RESIDUE, HEX_RESIDUE, NEUAC_RESIDUE, NEUGC_RESIDUE, DHEX_RESIDUE};

    // Common singly-charged oxonium B ions: {compact label, m/z}
    // m/z = sum(sugar_residue_masses) + PROTON
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

    // PPM tolerance for peak matching
    private static final double PPM_TOL = 20.0;

    /**
     * Annotate a spectrum.
     *
     * @param sequence       peptide sequence (single-letter codes, upper-case)
     * @param mods           list of modifications (position 1-indexed)
     * @param precursorCharge precursor charge state
     * @param specMzs        sorted observed m/z values
     * @param specInts       observed intensities (same order as specMzs)
     * @param ionTypes       list of ion type tokens: "a","b","c","x","y","z"
     * @param addGlycanIons  whether to generate glycan B/Y ions (for ori variants)
     * @param addNeutralLoss whether to annotate neutral loss variants
     * @return list of matched IonMatch objects
     */
    public static ArrayList<IonMatch> annotate(
            String sequence,
            ArrayList<ModificationMatch> mods,
            int precursorCharge,
            double[] specMzs,
            double[] specInts,
            List<String> ionTypes,
            boolean addGlycanIons,
            boolean addNeutralLoss) {

        ArrayList<IonMatch> result = new ArrayList<>();
        if (specMzs == null || specMzs.length == 0) return result;

        int n = sequence.length();

        // Build per-position mod mass arrays (1-indexed: index 1..n for residues)
        double nTermDelta = 0.0;
        double cTermDelta = 0.0;
        double[] resDelta = new double[n + 1]; // resDelta[1..n]

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
        double[] prefix = new double[n + 1]; // prefix[i] = sum of residue+mod masses for AA[0..i-1]
        prefix[0] = 0.0;
        for (int i = 1; i <= n; i++) {
            double aaMass = AA_MASS.getOrDefault(sequence.charAt(i - 1), 0.0);
            prefix[i] = prefix[i - 1] + aaMass + resDelta[i];
        }

        double totalPeptideMass = prefix[n] + H2O + nTermDelta + cTermDelta;

        // Max charge to generate (cap at precursorCharge, max 4)
        int maxCharge = Math.min(precursorCharge, 4);

        // ── Neutral losses to apply (when enabled) ──────────────────────────
        List<double[]> neutralLosses = new ArrayList<>(); // {mass, name_suffix}
        if (addNeutralLoss) {
            neutralLosses.add(new double[]{H2O, 0}); // use index trick; names built below
            neutralLosses.add(new double[]{NH3, 1});
            for (ModificationMatch mm : mods) {
                if (isPhospho(mm)) {
                    neutralLosses.add(new double[]{H3PO4, 2});
                    break; // add H3PO4 once
                }
            }
        }
        String[] lossNames = {"-H2O", "-NH3", "-H3PO4"};

        // ── Peptide fragment ions ────────────────────────────────────────────
        boolean doA = ionTypes.contains("a");
        boolean doB = ionTypes.contains("b");
        boolean doC = ionTypes.contains("c");
        boolean doX = ionTypes.contains("x");
        boolean doY = ionTypes.contains("y");
        boolean doZ = ionTypes.contains("z");

        for (int i = 1; i < n; i++) { // cleavage after position i (1-indexed)

            // b[i] neutral mass
            double bNeutral = prefix[i] + nTermDelta;
            // y[n-i] neutral mass
            double yNeutral = prefix[n] - prefix[i] + H2O + cTermDelta;

            for (int z = 1; z <= maxCharge; z++) {
                String zSuffix = chargeStr(z);

                if (doA) addMatch(result, specMzs, specInts,
                        "a" + i + zSuffix, (bNeutral - CO + z * PROTON) / z);
                if (doB) addMatch(result, specMzs, specInts,
                        "b" + i + zSuffix, (bNeutral + z * PROTON) / z);
                if (doC) addMatch(result, specMzs, specInts,
                        "c" + i + zSuffix, (bNeutral + NH3 + z * PROTON) / z);
                // x ion: y + CO2 - H2O
                if (doX) addMatch(result, specMzs, specInts,
                        "x" + (n - i) + zSuffix, (yNeutral - H2O + CO2 + z * PROTON) / z);
                if (doY) addMatch(result, specMzs, specInts,
                        "y" + (n - i) + zSuffix, (yNeutral + z * PROTON) / z);
                // z• (radical): y - NH2• = y - NH3 + H_atom
                if (doZ) addMatch(result, specMzs, specInts,
                        "z" + (n - i) + zSuffix, (yNeutral - Z_DOT_DELTA + z * PROTON) / z);

                // Neutral losses on b and y — charge suffix BEFORE neutral loss
                for (double[] loss : neutralLosses) {
                    int li = (int) loss[1];
                    String lName = lossNames[li];
                    double lMass = loss[0];
                    if (doB) addMatch(result, specMzs, specInts,
                            "b" + i + zSuffix + lName, (bNeutral - lMass + z * PROTON) / z);
                    if (doY) addMatch(result, specMzs, specInts,
                            "y" + (n - i) + zSuffix + lName, (yNeutral - lMass + z * PROTON) / z);
                }
            }
        }

        // ── Glycan ions (only for ori annotation passes) ────────────────────
        if (addGlycanIons) {
            // Identify glycan mods (N-glycan: mass>500 on N; O-glycan: mass>100 on S/T)
            double glycanMass = 0.0;
            boolean isNGlycan = false;
            for (ModificationMatch mm : mods) {
                String ptm = mm.getTheoreticPtm();
                if (ptm.contains(" of ")) {
                    String residue = ptm.split(" of ")[1];
                    double mass = mm.getMass();
                    if (residue.equals("N") && (mass > 500 || ptm.startsWith("203.079") || ptm.startsWith("0.0"))) {
                        glycanMass = mass;
                        isNGlycan = true;
                        break;
                    } else if ((residue.equals("S") || residue.equals("T") || residue.equals("ST"))
                            && (mass > 100 || ptm.startsWith("0.0"))) {
                        glycanMass = mass;
                        break;
                    }
                }
            }

            // Glycan B (oxonium) ions — singly charged only
            for (Object[] entry : OXONIUM_IONS) {
                double mz = (double) entry[1];
                int peakIdx = findBestPeak(specMzs, mz);
                if (peakIdx >= 0)
                    result.add(new IonMatch((String) entry[0], specMzs[peakIdx], specInts[peakIdx], mz));
            }

            // Glycan Y ions: bare peptide + partial glycan compositions
            // Bare peptide mass = total mass - glycan mass at N/O site
            if (glycanMass > 0) {
                double barePeptideMass = totalPeptideMass - glycanMass;
                generateGlycanYIons(result, specMzs, specInts, barePeptideMass, glycanMass, maxCharge);
            }
        }

        return result;
    }

    /**
     * Generate glycan Y ions by enumerating partial glycan compositions
     * (as sugar residue combinations) up to the full glycan mass.
     */
    private static void generateGlycanYIons(
            ArrayList<IonMatch> result,
            double[] specMzs, double[] specInts,
            double barePeptideMass,
            double glycanMass,
            int maxCharge) {

        // Max amounts of each sugar (HexNAc, Hex, NeuAc, NeuGc, dHex)
        // Capped at reasonable values; the mass guard will prune over-large combos.
        int[] maxCounts = {10, 10, 5, 3, 4};

        // Depth-limited enumeration — 5 nested loops but pruned early by mass.
        // Only generate combos with total residue mass <= glycanMass + 1 Da.
        // Skip the empty combo (0,0,0,0,0).
        for (int n0 = 0; n0 <= maxCounts[0]; n0++) {
            double m0 = n0 * SUGAR_MASSES[0];
            if (m0 > glycanMass + 1.0) break;
            for (int n1 = 0; n1 <= maxCounts[1]; n1++) {
                double m1 = m0 + n1 * SUGAR_MASSES[1];
                if (m1 > glycanMass + 1.0) break;
                for (int n2 = 0; n2 <= maxCounts[2]; n2++) {
                    double m2 = m1 + n2 * SUGAR_MASSES[2];
                    if (m2 > glycanMass + 1.0) break;
                    for (int n3 = 0; n3 <= maxCounts[3]; n3++) {
                        double m3 = m2 + n3 * SUGAR_MASSES[3];
                        if (m3 > glycanMass + 1.0) break;
                        for (int n4 = 0; n4 <= maxCounts[4]; n4++) {
                            double comboMass = m3 + n4 * SUGAR_MASSES[4];
                            if (comboMass > glycanMass + 1.0) break;
                            if (n0 + n1 + n2 + n3 + n4 == 0) continue; // skip bare peptide

                            double yMass = barePeptideMass + comboMass;
                            // Build label lazily — only when a peak match is found.
                            String label = null;
                            for (int z = 1; z <= maxCharge; z++) {
                                double mz = (yMass + z * PROTON) / z;
                                int peakIdx = findBestPeak(specMzs, mz);
                                if (peakIdx >= 0) {
                                    if (label == null) label = buildGlycanYLabel(n0, n1, n2, n3, n4);
                                    result.add(new IonMatch(label + chargeStr(z), specMzs[peakIdx], specInts[peakIdx], mz));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Build a compact Y_ label from sugar counts.
     * e.g. (2, 1, 0, 0, 0) → "Y_N2H1"
     */
    private static String buildGlycanYLabel(int... counts) {
        StringBuilder sb = new StringBuilder("Y_");
        for (int i = 0; i < SUGAR_LETTERS.length; i++) {
            if (counts[i] > 0) sb.append(SUGAR_LETTERS[i]).append(counts[i]);
        }
        return sb.toString();
    }

    // Precomputed charge suffix strings to avoid repeated String.repeat() allocations
    private static final String[] CHARGE_STR = {"", "", "++", "+++", "++++"};


    /** Try to match a theoretical m/z to the nearest observed peak within PPM_TOL. */
    private static void addMatch(
            ArrayList<IonMatch> result,
            double[] specMzs, double[] specInts,
            String label, double theoMz) {

        if (theoMz <= 0) return;
        int best = findBestPeak(specMzs, theoMz);
        if (best >= 0) result.add(new IonMatch(label, specMzs[best], specInts[best], theoMz));
    }

    /**
     * Returns the index of the nearest peak within PPM_TOL, or -1 if none.
     * Used by the glycan Y-ion loop to defer label building until a hit is confirmed.
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
        double ppm = Math.abs(specMzs[best] - theoMz) / theoMz * 1e6;
        return ppm <= PPM_TOL ? best : -1;
    }

    /** Returns charge indicator string: "" for z=1, "++" for z=2, etc. */
    private static String chargeStr(int z) {
        return z < CHARGE_STR.length ? CHARGE_STR[z] : "+".repeat(z);
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
