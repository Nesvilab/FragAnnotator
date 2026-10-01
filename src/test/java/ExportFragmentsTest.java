import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

/** Pins that the appended ion columns always start at the same column as their headers. */
class ExportFragmentsTest {

    private static final String HEADER =
            "Spectrum\tPeptide\tGene\tProtein Description\tMapped Genes\tMapped Proteins";

    /** The column the first ion value lands in when this row is written the way the tool writes it. */
    private static int firstIonColumn(String[] aligned) {
        return Arrays.asList((String.join("\t", aligned) + "\tb2").split("\t", -1)).indexOf("b2");
    }

    /** The column the first ion header lands in. */
    private static int firstIonHeaderColumn(String headerLine) {
        return Arrays.asList((headerLine.stripTrailing() + "\tions").split("\t", -1)).indexOf("ions");
    }

    @Test
    void aRowMissingItsTrailingProteinFieldsIsPaddedToTheHeader() {
        // A writer that left off Gene..Mapped Proteins instead of padding them with empty fields.
        String[] row = "s.1.1.2\tPEPTIDE".split("\t", -1);
        int width = ExportFragments.headerWidth(HEADER);
        String[] aligned = ExportFragments.alignToHeader(row, width);
        assertEquals(6, aligned.length);
        assertEquals("PEPTIDE", aligned[1]);
        assertEquals("", aligned[5]);
        assertEquals(firstIonHeaderColumn(HEADER), firstIonColumn(aligned));
    }

    @Test
    void aFullWidthRowIsUnchanged() {
        String[] row = "s.1.1.2\tPEPTIDE\t\t\t\t".split("\t", -1);
        String[] aligned = ExportFragments.alignToHeader(row, ExportFragments.headerWidth(HEADER));
        assertSame(row, aligned);
        assertEquals(firstIonHeaderColumn(HEADER), firstIonColumn(aligned));
    }

    @Test
    void aTrailingTabOnTheHeaderIsNotAColumn() {
        // The header is written back without trailing whitespace, so counting the empty field after
        // a trailing tab would put every row's ion columns one place right of their headers.
        String header = HEADER + "\t";
        int width = ExportFragments.headerWidth(header);
        assertEquals(6, width);
        String[] aligned = ExportFragments.alignToHeader("s.1.1.2\tPEPTIDE".split("\t", -1), width);
        assertEquals(firstIonHeaderColumn(header), firstIonColumn(aligned));
    }

    @Test
    void emptyFieldsPastTheHeaderAreDropped() {
        // A writer that ends every row with a tab the header does not have.
        String[] row = "s.1.1.2\tPEPTIDE\tG\tD\tMG\tMP\t\t".split("\t", -1);
        String[] aligned = ExportFragments.alignToHeader(row, ExportFragments.headerWidth(HEADER));
        assertEquals(6, aligned.length);
        assertEquals("MP", aligned[5]);
        assertEquals(firstIonHeaderColumn(HEADER), firstIonColumn(aligned));
    }

    @Test
    void anEmptyInteriorHeaderStillCounts() {
        String header = "Spectrum\t\tPeptide";
        assertEquals(3, ExportFragments.headerWidth(header));
    }

    @Test
    void aValuePastTheLastHeaderCannotBeAligned() {
        // Placing it would mean discarding data or putting an ion column under the wrong header.
        String[] row = "s.1.1.2\tPEPTIDE\tG\tD\tMG\tMP\textra".split("\t", -1);
        assertNull(ExportFragments.alignToHeader(row, ExportFragments.headerWidth(HEADER)));
    }
}
