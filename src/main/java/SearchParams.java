import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.*;

/**
 * The MSFragger search settings read out of {@code fragpipe.workflow}.
 *
 * <p><b>This is the one reader of that file.</b> A FragPipe result folder records the search it
 * came from, so the annotator uses the ion series, fragment tolerance, custom series and labile
 * mass offsets that search actually used instead of its own defaults. Defaults of b/y at 20 ppm
 * are wrong for most real searches — on an EThcD or a labile run they are precisely the ions that
 * did not identify the peptide.
 *
 * <p>Nothing here is inferred from the data. A parameter that is absent leaves its field at the
 * default, because annotating ion series the search never generated is worse than annotating none.
 */
public class SearchParams {

    /** {@code msfragger.fragment_ion_series} verbatim, case preserved (uppercase Y is not y). */
    public List<String> ionSeries = new ArrayList<>();
    /** Fragment tolerance; {@link #DEFAULT_TOL_PPM} ppm when the parameter is absent. */
    public double fragTol = DEFAULT_TOL_PPM;
    /** True when {@link #fragTol} is in Da rather than ppm ({@code fragment_mass_units=0}). */
    public boolean fragTolDa = false;
    /** User-defined series from {@code msfragger.ion_series_definitions}. */
    public List<CustomIon> customIons = new ArrayList<>();
    /** {@code msfragger.labile_fragment_ion_series} — the base series remainders were generated for. */
    public List<String> labileSeries = new ArrayList<>();
    /** Whether the search ran in a labile mode at all. Decides whether the labile columns exist. */
    public boolean labile = false;
    /** True when {@code labile_search_mode=nglycan}; used only to warn about a missing glycan column. */
    public boolean nglycanMode = false;
    /** The declared labile mass offsets. */
    public List<LabileOffset> offsets = new ArrayList<>();
    /**
     * The search's global labile ion lists as a single site-less offset, or null when it declared
     * none. Used for a modification no declared offset explains — above all a glycan, because a
     * glyco search keeps its glycan masses in {@code mass_offsets_detailed} and leaves
     * {@code mass_offsets=0}, so with detailed offsets switched off there is no per-mass entry to
     * find and the glycan would otherwise get no labile ions at all.
     */
    public LabileOffset globalOffset = null;

    /** Fallback tolerance for a result whose workflow does not state one (the previous hardcoded value). */
    public static final double DEFAULT_TOL_PPM = 20.0;

    /**
     * How far a psm.tsv modification mass may sit from a declared offset and still be it. psm.tsv
     * writes 4 decimals, so rounding alone is under 5e-5; this window covers that while staying
     * narrow enough not to merge two neighbouring offsets.
     */
    public static final double OFFSET_TOL = 0.005;

    // Reading ---------------------------------------------------------------

    /** Parses {@code <dir>/fragpipe.workflow}. A missing file yields defaults, which is normal. */
    public static SearchParams read(File dir) {
        File f = new File(dir, "fragpipe.workflow");
        if (!f.isFile()) {
            System.out.println("No fragpipe.workflow found; using default ion types and "
                    + DEFAULT_TOL_PPM + " ppm tolerance.");
            return new SearchParams();
        }
        try {
            return parse(new String(Files.readAllBytes(f.toPath()), StandardCharsets.UTF_8));
        } catch (IOException e) {
            System.err.println("WARNING: could not read " + f + ": " + e.getMessage());
            return new SearchParams();
        }
    }

    /** Parses workflow-file text. Package-visible so tests can pin the format without a file. */
    static SearchParams parse(String text) {
        Map<String, String> p = new HashMap<>();
        for (String raw : text.split("\r?\n")) {
            int lead = 0;
            while (lead < raw.length() && Character.isWhitespace(raw.charAt(lead))) lead++;
            String line = raw.substring(lead);
            if (line.isEmpty() || line.charAt(0) == '#' || line.charAt(0) == '!') continue;
            // A key never contains an escaped '=', so the first '=' is always the separator.
            int eq = line.indexOf('=');
            if (eq < 0) continue;
            p.put(line.substring(0, eq).trim(), unescape(line.substring(eq + 1)));
        }

        SearchParams s = new SearchParams();
        s.ionSeries = tokens(get(p, "msfragger.fragment_ion_series"));
        s.customIons = parseCustomIons(get(p, "msfragger.ion_series_definitions"));
        s.labileSeries = tokens(get(p, "msfragger.labile_fragment_ion_series"));

        Double tol = parseDouble(get(p, "msfragger.fragment_mass_tolerance"));
        if (tol != null) {
            s.fragTol = tol;
            // 1 = ppm, 0 = Da. Anything else (including an absent parameter) stays ppm, which is
            // both the MSFragger default and this tool's previous behaviour.
            String units = get(p, "msfragger.fragment_mass_units");
            s.fragTolDa = units.equals("0") || units.equalsIgnoreCase("da");
        }

        String mode = get(p, "msfragger.labile_search_mode");
        s.labile = !(mode.isEmpty() || mode.equalsIgnoreCase("off"));
        s.nglycanMode = mode.equalsIgnoreCase("nglycan");

        if (s.labile) {
            // The detailed list survives in the workflow of a search that switched labile mode off
            // or never enabled detailed offsets, and those offsets were never searched.
            if (truthy(get(p, "msfragger.use_detailed_offsets"))) {
                String raw = p.containsKey("msfragger.mass_offsets_detailed")
                        ? p.get("msfragger.mass_offsets_detailed")
                        : get(p, "msfragger.detailed_mass_offsets");
                for (String entry : raw.split(";")) {
                    LabileOffset off = LabileOffset.parseEntry(entry);
                    if (off != null) s.offsets.add(off);
                }
            } else {
                s.offsets = globalOffsets(p);
            }
            s.globalOffset = globalOffset(p);
        }
        return s;
    }

    /**
     * The labile offsets of a search that did not use the detailed list: MSFragger then applies
     * three global ion lists to every declared mass offset. Without this a labile search using the
     * simple parameters would produce three empty columns — present, so the file looks annotated.
     */
    private static List<LabileOffset> globalOffsets(Map<String, String> p) {
        LabileOffset global = globalOffset(p);
        List<LabileOffset> out = new ArrayList<>();
        if (global == null) return out;
        for (double mass : numbers(get(p, "msfragger.mass_offsets"))) {
            if (Math.abs(mass) < 1e-6) continue; // the mandatory no-offset entry
            out.add(new LabileOffset(mass, "", global.diagnostic, global.peptide, global.fragment));
        }
        return out;
    }

    /**
     * The three global ion lists as one site-less offset, or null when the search declares none.
     * Its mass is meaningless and it is never matched by mass — it is the ion lists for a
     * modification that has no entry of its own.
     */
    private static LabileOffset globalOffset(Map<String, String> p) {
        double[] diag = numbers(get(p, "msfragger.diagnostic_fragments"));
        double[] pep = numbers(get(p, "msfragger.Y_type_masses"));
        double[] frag = numbers(get(p, "msfragger.remainder_fragment_masses"));
        if (diag.length == 0 && pep.length == 0 && frag.length == 0) return null;
        return new LabileOffset(0.0, "", diag, pep, frag);
    }

    // Lookup ----------------------------------------------------------------

    /**
     * The declared offset explaining a modification of {@code mass} on residue {@code aa}, or null.
     *
     * <p>The site matters, not just the mass: MSFragger's list routinely declares the same mass
     * twice with different sites and different remainder ions. The ADPr reference search has
     * 541.0611 on {@code SKTYHDE} (fragment remainder 0) and on {@code R} (-42.0205), so matching
     * on mass alone annotates the wrong ions for the residue at hand.
     *
     * @param aa the modified residue, or 0 when the modification is unlocalized and has none
     */
    public LabileOffset offsetFor(double mass, char aa) {
        if (Math.abs(mass) < 1e-6) return null;
        LabileOffset best = null;
        double bestDiff = Double.MAX_VALUE;
        for (LabileOffset o : offsets) {
            double diff = Math.abs(o.mass - mass);
            if (diff > OFFSET_TOL) continue;
            // An unlocalized modification has no residue to test, so the site check cannot run.
            if (aa != 0 && !o.allows(aa)) continue;
            if (diff < bestDiff) {
                bestDiff = diff;
                best = o;
            }
        }
        return best;
    }

    /** True when any declared offset states diagnostic ions. */
    public boolean hasDiagnostic() {
        for (LabileOffset o : offsets) if (o.diagnostic.length > 0) return true;
        return false;
    }

    /** True when any declared offset states peptide remainders. */
    public boolean hasPepRemainder() {
        for (LabileOffset o : offsets) if (o.peptide.length > 0) return true;
        return false;
    }

    // Text helpers ----------------------------------------------------------

    private static String get(Map<String, String> p, String k) {
        String v = p.get(k);
        return v == null ? "" : v.trim();
    }

    /**
     * Java {@code .properties} unescaping: a backslash escapes the character after it, with the
     * usual control escapes. Without this every number after a {@code \=} in
     * {@code mass_offsets_detailed} keeps its backslash and parses as nothing at all.
     */
    static String unescape(String s) {
        StringBuilder out = new StringBuilder(s.length());
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            if (c != '\\') {
                out.append(c);
                continue;
            }
            if (++i >= s.length()) break;
            char n = s.charAt(i);
            switch (n) {
                case 'n': out.append('\n'); break;
                case 't': out.append('\t'); break;
                case 'r': out.append('\r'); break;
                default:  out.append(n);
            }
        }
        return out.toString();
    }

    /**
     * Splits a list-valued parameter, tolerating every separator FragPipe has written one with:
     * commas, semicolons, slashes and whitespace. Both {@code b,y} and
     * {@code 0/114.03169/193.99802} are real, and {@code ion_series_definitions} has been written
     * both {@code ;}- and space-separated.
     */
    static List<String> tokens(String v) {
        List<String> out = new ArrayList<>();
        if (v == null) return out;
        for (String t : v.split("[,;/\\s]+")) {
            t = t.trim();
            if (!t.isEmpty()) out.add(t);
        }
        return out;
    }

    /** {@link #tokens} parsed as numbers, silently dropping anything that is not one. */
    static double[] numbers(String v) {
        List<String> toks = tokens(v);
        double[] tmp = new double[toks.size()];
        int n = 0;
        for (String t : toks) {
            Double d = parseDouble(t);
            if (d != null) tmp[n++] = d;
        }
        return Arrays.copyOf(tmp, n);
    }

    static Double parseDouble(String s) {
        if (s == null || s.trim().isEmpty()) return null;
        try {
            return Double.parseDouble(s.trim());
        } catch (NumberFormatException e) {
            return null;
        }
    }

    /** FragPipe writes booleans as true/false; MSFragger's own params file uses 1/0. */
    static boolean truthy(String v) {
        return v.equalsIgnoreCase("true") || v.equals("1");
    }

    /**
     * {@code msfragger.ion_series_definitions} is a flat run of (name, terminus, offset) triples,
     * whatever it is separated by — FragPipe has written the same definitions {@code ;}-separated
     * in the workflow and space-separated in its own log.
     */
    static List<CustomIon> parseCustomIons(String v) {
        List<String> t = tokens(v);
        List<CustomIon> out = new ArrayList<>();
        for (int i = 0; i + 2 < t.size(); i += 3) {
            String name = t.get(i);
            String term = t.get(i + 1);
            Double off = parseDouble(t.get(i + 2));
            if (name.isEmpty() || off == null) continue;
            boolean nterm;
            if (term.equalsIgnoreCase("N")) nterm = true;
            else if (term.equalsIgnoreCase("C")) nterm = false;
            else continue;
            out.add(new CustomIon(name, nterm, off));
        }
        return out;
    }
}
