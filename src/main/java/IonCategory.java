/**
 * The four kinds of ion this tool annotates. Each becomes one group of five columns in psm.tsv.
 *
 * <p>Backbone ions are always annotated. The other three exist only for a labile search, and are
 * the generic form of what the glyco-specific annotation used to produce: diagnostic ions are the
 * oxonium markers, peptide remainders the Y ions, and fragment remainders the backbone remnants.
 */
public enum IonCategory {

    /** Standard and custom series with all modifications intact. */
    BACKBONE(""),
    /** Backbone ions whose labile modification was partly or wholly lost. Positional. */
    FRAG_REMAINDER("frag_rem_"),
    /** The intact peptide carrying a partial modification loss. Not positional. */
    PEP_REMAINDER("pep_rem_"),
    /** Modification-specific marker ions detached from the peptide. Not positional. */
    DIAGNOSTIC("diag_");

    /** Column-name prefix; empty for the backbone group, which keeps the original column names. */
    public final String prefix;

    IonCategory(String prefix) {
        this.prefix = prefix;
    }
}
