import java.io.*;
import java.util.*;

public class FragmentsExportMainClass {

    /** Recognised backbone ion letters. */
    private static final Set<String> LETTERS =
            new LinkedHashSet<>(Arrays.asList("a", "b", "c", "x", "y", "z"));
    /**
     * Ion-type keywords the glyco modes used. They no longer select anything: glycan ions are now
     * annotated from the search's own labile parameters, and whether a PSM gets composition labels
     * is decided by its own glycan composition. Accepted and ignored so an older FragPipe that
     * still passes them does not fail.
     */
    private static final Set<String> RETIRED_KEYWORDS =
            new LinkedHashSet<>(Arrays.asList("ngly", "ogly", "gly"));

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
        Options options = parseIonTypes(ionsTypes);
        new ExportFragments(resultsFolder, threadsNumber, options.letters, options.neutralLoss,
                glycanResiduesPath, glycanModsPath);
    }

    /** What the ion-types argument asked for, beyond what the search itself states. */
    static class Options {
        /** Backbone letters to annotate instead of the search's own; empty means "use the search's". */
        final List<String> letters;
        /** Whether to annotate neutral-loss variants of b and y. */
        final boolean neutralLoss;

        Options(List<String> letters, boolean neutralLoss) {
            this.letters = letters;
            this.neutralLoss = neutralLoss;
        }
    }

    /**
     * Parse the ion-types argument.
     *
     * <p>The search's own {@code fragpipe.workflow} is the source of truth for which ion series
     * were generated, at what tolerance, and with what custom and labile definitions. This
     * argument only overrides the standard backbone letters, for a user who wants a different
     * picture from the one the search produced.
     *
     * <p>Contract with FragPipe (see {@code CmdExportMatchedFragments.java} in the FragPipe repo):
     * when the user leaves the fragment-type selector at its default, FragPipe passes the sentinel
     * {@code "r"}. That — like an empty argument — means "take the ion series from the search",
     * which is now the case that does the right thing on its own.
     *
     * <p>Accepts a flexible input format so hand-typed args are not silently punished:
     * <ul>
     *   <li>separators: {@code _}, {@code ,}, whitespace, {@code +}, {@code /}</li>
     *   <li>backbone ion letters: {@code a b c x y z}, also concatenated ({@code by})</li>
     *   <li>{@code neu}: also annotate neutral losses of b and y</li>
     *   <li>{@code ngly}/{@code ogly}/{@code gly}: accepted and ignored (see
     *       {@link #RETIRED_KEYWORDS})</li>
     * </ul>
     */
    static Options parseIonTypes(String ionsTypes) {
        LinkedHashSet<String> letters = new LinkedHashSet<>();
        boolean neutralLoss = false;
        List<String> unrecognized = new ArrayList<>();
        List<String> retired = new ArrayList<>();

        if (ionsTypes != null && ionsTypes.trim().equalsIgnoreCase("r")) {
            System.out.println("Ion types: taking the search's own fragment_ion_series "
                    + "(FragPipe default selection).");
            return new Options(Collections.emptyList(), false);
        }

        if (ionsTypes != null && !ionsTypes.trim().isEmpty()) {
            for (String token : ionsTypes.trim().toLowerCase().split("[_,\\s+/]+")) {
                if (token.isEmpty()) continue;
                if (token.equals("neu")) {
                    neutralLoss = true;
                } else if (RETIRED_KEYWORDS.contains(token)) {
                    retired.add(token);
                } else if (LETTERS.contains(token)) {
                    letters.add(token);
                } else if (allLetters(token)) {
                    // Concatenated backbone letters (e.g. "by" -> b,y)
                    for (int i = 0; i < token.length(); i++) letters.add(String.valueOf(token.charAt(i)));
                } else {
                    unrecognized.add(token);
                }
            }
        }

        if (!retired.isEmpty()) {
            System.out.println("Note: ion-type keyword(s) " + retired + " are no longer needed and "
                    + "were ignored. Glycan and other labile ions are now annotated automatically "
                    + "from the search's own parameters.");
        }
        if (!unrecognized.isEmpty()) {
            System.err.println("WARNING: unrecognized ion-type token(s): " + unrecognized
                    + ". Valid letters: a/b/c/x/y/z. Valid keyword: neu.");
        }

        List<String> list = new ArrayList<>(letters);
        if (list.isEmpty()) {
            System.out.println("Ion types: taking the search's own fragment_ion_series.");
        } else {
            System.out.println("Ion types: " + list + " (overriding the search's fragment_ion_series)");
        }
        return new Options(list, neutralLoss);
    }

    private static boolean allLetters(String token) {
        for (int i = 0; i < token.length(); i++) {
            if (!LETTERS.contains(String.valueOf(token.charAt(i)))) return false;
        }
        return !token.isEmpty();
    }
}
