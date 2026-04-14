import java.io.*;
import java.util.*;

public class FragmentsExportMainClass {

    // Recognised multi-character keywords (non-backbone modifiers).
    private static final Set<String> KEYWORDS =
            new LinkedHashSet<>(Arrays.asList("ngly", "ogly", "gly", "neu"));
    // Recognised backbone ion letters.
    private static final Set<String> LETTERS  =
            new LinkedHashSet<>(Arrays.asList("a", "b", "c", "x", "y", "z"));

    public static void main(String[] args) {
        File resultsFolder    = new File(args[0]);
        int threadsNumber     = Integer.parseInt(args[1]);
        String ionsTypes      = args[2];
        String glycanResiduesPath = args.length > 3 ? args[3] : "";
        String glycanModsPath     = args.length > 4 ? args[4] : "";
        try {
            new FragmentsExportMainClass(resultsFolder, threadsNumber, ionsTypes,
                    glycanResiduesPath, glycanModsPath);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public FragmentsExportMainClass(File resultsFolder, int threadsNumber, String ionsTypes,
            String glycanResiduesPath, String glycanModsPath) throws IOException {
        List<String> ionsTypeArray = parseIonTypes(ionsTypes);
        new ExportFragments(resultsFolder, threadsNumber, ionsTypeArray,
                glycanResiduesPath, glycanModsPath);
    }

    /**
     * Parse the ion-types argument into a canonical list of tokens.
     *
     * <p>Contract with FragPipe (see {@code CmdExportMatchedFragments.java}
     * in the FragPipe repo): when the user leaves the fragment-type selector
     * at its default (b and y), FragPipe passes the sentinel string
     * {@code "r"} on the command line. This parser therefore treats
     * {@code "r"} (case-insensitive) as an alias for {@code b_y}.
     *
     * <p>Additionally accepts a flexible input format so hand-typed args are
     * not silently punished for using commas, spaces, or concatenated letters:
     * <ul>
     *   <li>separators: {@code _}, {@code ,}, whitespace, {@code +}, {@code /}</li>
     *   <li>multi-char keywords (case-insensitive): {@code ngly}, {@code ogly},
     *       {@code gly}, {@code neu}</li>
     *   <li>backbone ion letters: {@code a b c x y z}</li>
     *   <li>any remaining token whose characters are all backbone letters is
     *       split per-character, so {@code by} &rarr; {@code [b, y]}</li>
     * </ul>
     *
     * <p>If after parsing there is no backbone ion letter (a/b/c/x/y/z), a
     * warning is printed and the parser defaults to {@code b, y} so the run
     * still produces useful annotations instead of silently emitting empty
     * columns.
     */
    static List<String> parseIonTypes(String ionsTypes) {
        LinkedHashSet<String> out = new LinkedHashSet<>();
        List<String> unrecognized = new ArrayList<>();

        // FragPipe default-selection sentinel → b, y
        if (ionsTypes != null && ionsTypes.trim().equalsIgnoreCase("r")) {
            out.add("b");
            out.add("y");
            List<String> list = new ArrayList<>(out);
            System.out.println("Using ion types: " + list + " (FragPipe default)");
            return list;
        }

        if (ionsTypes != null && !ionsTypes.trim().isEmpty()) {
            String[] tokens = ionsTypes.trim().toLowerCase().split("[_,\\s+/]+");
            for (String token : tokens) {
                if (token.isEmpty()) continue;
                if (KEYWORDS.contains(token) || LETTERS.contains(token)) {
                    out.add(token);
                    continue;
                }
                // Concatenated backbone letters (e.g. "by" -> b,y)
                boolean allLetters = !token.isEmpty();
                for (int i = 0; i < token.length(); i++) {
                    if (!LETTERS.contains(String.valueOf(token.charAt(i)))) {
                        allLetters = false;
                        break;
                    }
                }
                if (allLetters) {
                    for (int i = 0; i < token.length(); i++) {
                        out.add(String.valueOf(token.charAt(i)));
                    }
                } else {
                    unrecognized.add(token);
                }
            }
        }

        if (!unrecognized.isEmpty()) {
            System.err.println("WARNING: unrecognized ion-type token(s): " + unrecognized
                    + ". Valid letters: a/b/c/x/y/z. Valid keywords: ngly/ogly/gly/neu.");
        }

        // Ensure at least one backbone ion letter is present, else default to b,y.
        boolean hasBackbone = false;
        for (String t : out) {
            if (LETTERS.contains(t)) { hasBackbone = true; break; }
        }
        if (!hasBackbone) {
            System.err.println("WARNING: no valid backbone ion types parsed from '" + ionsTypes
                    + "'; defaulting to b, y");
            out.add("b");
            out.add("y");
        }

        List<String> list = new ArrayList<>(out);
        System.out.println("Using ion types: " + list);
        return list;
    }
}
