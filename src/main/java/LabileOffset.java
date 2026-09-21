import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * One entry of {@code msfragger.mass_offsets_detailed} — a labile modification mass, the residues
 * it is allowed on, and the three kinds of ion it can produce:
 *
 * <ul>
 *   <li>{@code _d=} <b>diagnostic</b> — low-m/z marker ions, given as <b>m/z</b> (the glyco
 *       oxonium ions' analogue).</li>
 *   <li>{@code _p=} <b>peptide remainder</b> — masses the modification leaves on the intact
 *       peptide (the glycan Y ions' analogue).</li>
 *   <li>{@code _f=} <b>fragment remainder</b> — masses it leaves on a backbone fragment that
 *       contains it.</li>
 * </ul>
 *
 * <p>A remainder of 0 means the modification was lost entirely, leaving the bare residue. It is
 * added to both remainder lists whether or not the search declared it, because complete loss is
 * physically available to any labile modification and the declared lists routinely omit it — a
 * glyco offset declares only {@code _f=203.07937}, so without this the backbone ions of a
 * glycopeptide with the glycan fully stripped would never be annotated.
 */
public class LabileOffset {

    public final double mass;
    /** The {@code aa=} value verbatim, kept for diagnostics. */
    public final String rawSites;
    /** Diagnostic ion m/z values, exactly as declared. */
    public final double[] diagnostic;
    /** Peptide remainder masses, including 0. */
    public final double[] peptide;
    /** Fragment remainder masses, including 0. */
    public final double[] fragment;

    /** The residues this offset may sit on. Empty means any residue. */
    private final char[] sites;
    /** True when {@link #sites} is a negated class ({@code [^P]}): any residue EXCEPT those. */
    private final boolean negated;

    public LabileOffset(double mass, String rawSites, double[] diagnostic,
                        double[] peptide, double[] fragment) {
        this.mass = mass;
        this.rawSites = rawSites == null ? "" : rawSites;
        this.diagnostic = dedup(diagnostic, false);
        this.peptide = dedup(peptide, true);
        this.fragment = dedup(fragment, true);

        char[] parsed = parseSites(this.rawSites);
        this.negated = parsed.length > 0 && parsed[0] == '^';
        this.sites = negated ? Arrays.copyOfRange(parsed, 1, parsed.length) : parsed;
    }

    /** Whether this offset is allowed on residue {@code aa}. An empty site list means any. */
    public boolean allows(char aa) {
        if (sites.length == 0) return true;
        boolean found = false;
        char up = Character.toUpperCase(aa);
        for (char s : sites) {
            if (Character.toUpperCase(s) == up) {
                found = true;
                break;
            }
        }
        return negated != found;
    }

    /** True when this offset declares no ions at all, and so is not labile. */
    public boolean isEmpty() {
        return diagnostic.length == 0 && peptide.length == 0 && fragment.length == 0;
    }

    // Parsing ---------------------------------------------------------------

    /**
     * Parses one {@code mass(aa=..._d=..._p=..._f=...)} entry. Returns null for an entry that
     * declares no ions at all — the mandatory {@code 0.0000(aa=)} entry, and any offset with no
     * labile fragmentation.
     */
    public static LabileOffset parseEntry(String entry) {
        if (entry == null) return null;
        entry = entry.trim();
        if (entry.isEmpty()) return null;

        String massPart = entry;
        String body = "";
        int open = entry.indexOf('(');
        if (open >= 0) {
            // Close on the LAST ')', not the first: a site may itself contain parentheses
            // (a sequon naming its modified residue), which would otherwise truncate the body.
            int close = entry.lastIndexOf(')');
            massPart = entry.substring(0, open);
            body = close > open ? entry.substring(open + 1, close) : entry.substring(open + 1);
        }
        Double mass = SearchParams.parseDouble(massPart);
        if (mass == null || Math.abs(mass) < 1e-6) return null;

        String sites = "";
        double[] d = new double[0], p = new double[0], f = new double[0];
        for (String part : body.split("_")) {
            int eq = part.indexOf('=');
            if (eq < 0) continue;
            String tag = part.substring(0, eq).trim();
            String val = part.substring(eq + 1).trim();
            switch (tag) {
                case "aa": sites = val; break;
                case "d":  d = SearchParams.numbers(val); break;
                case "p":  p = SearchParams.numbers(val); break;
                case "f":  f = SearchParams.numbers(val); break;
                default:   break;
            }
        }
        // Build with the RAW lists so an offset declaring nothing stays empty. The implicit zero
        // must not make an inert offset look labile.
        if (d.length == 0 && p.length == 0 && f.length == 0) return null;
        return new LabileOffset(mass, sites, d, p, f);
    }

    /**
     * The residues an {@code aa=} value designates as the modification site.
     *
     * <p>Two syntaxes occur. A plain residue list ({@code SKTYHDE}) is itself the site set. A
     * sequon ({@code {N[^P][ST]}}) describes the sequence context the offset was searched in, and
     * only its <b>modified residue</b> is a site: the residue in parentheses when the sequon names
     * one, and the first residue otherwise. The rest of the motif was the search's own constraint
     * and is already satisfied by the search having assigned the modification there, so re-testing
     * it could only reject what the search accepted — and would fail outright on a semi-tryptic
     * peptide whose sequon is completed by protein residues the peptide does not contain.
     *
     * @return the allowed residues, or a set whose first element is {@code '^'} for a negated
     *         class; empty for "any residue"
     */
    static char[] parseSites(String raw) {
        if (raw == null) return new char[0];
        String s = raw.trim();
        if (s.isEmpty() || s.contains("*")) return new char[0];

        int lp = s.indexOf('(');
        if (lp >= 0) {
            int rp = s.indexOf(')', lp);
            if (rp > lp) return classChars(s.substring(lp + 1, rp));
        }
        if (s.startsWith("{")) {
            String motif = s.substring(1, s.endsWith("}") ? s.length() - 1 : s.length());
            return classChars(firstToken(motif));
        }
        return classChars(s);
    }

    /** The first motif token: a bracketed character class, or a single character. */
    private static String firstToken(String motif) {
        if (motif.isEmpty()) return "";
        if (motif.charAt(0) == '[') {
            int close = motif.indexOf(']');
            return close > 0 ? motif.substring(0, close + 1) : motif;
        }
        return motif.substring(0, 1);
    }

    /**
     * The residues a token designates. {@code [ST]} is S and T; {@code [^P]} is anything but P,
     * returned with a leading {@code '^'} marker; a bare run of letters is itself the set.
     */
    private static char[] classChars(String token) {
        String t = token.trim();
        boolean neg = false;
        if (t.startsWith("[") && t.endsWith("]") && t.length() >= 2) {
            t = t.substring(1, t.length() - 1);
        }
        if (t.startsWith("^")) {
            neg = true;
            t = t.substring(1);
        }
        List<Character> out = new ArrayList<>();
        if (neg) out.add('^');
        for (char c : t.toCharArray()) {
            if (Character.isLetter(c)) out.add(c);
        }
        // A negated class with nothing to exclude allows everything, which is the empty site set.
        if (neg && out.size() == 1) return new char[0];
        char[] arr = new char[out.size()];
        for (int i = 0; i < arr.length; i++) arr[i] = out.get(i);
        return arr;
    }

    /** Sorted unique masses, optionally with 0 added. */
    private static double[] dedup(double[] in, boolean withZero) {
        double[] tmp = new double[(in == null ? 0 : in.length) + (withZero ? 1 : 0)];
        int n = 0;
        if (withZero) tmp[n++] = 0.0;
        if (in != null) {
            for (double v : in) {
                boolean seen = false;
                for (int i = 0; i < n; i++) {
                    if (Math.abs(tmp[i] - v) < 1e-6) {
                        seen = true;
                        break;
                    }
                }
                if (!seen) tmp[n++] = v;
            }
        }
        double[] out = Arrays.copyOf(tmp, n);
        Arrays.sort(out);
        return out;
    }

    @Override
    public String toString() {
        return mass + "(aa=" + rawSites + " d=" + diagnostic.length
                + " p=" + peptide.length + " f=" + fragment.length + ")";
    }
}
