import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/** Pins the command-line contract with FragPipe. */
class FragmentsExportMainClassTest {

    @Test
    void theFragPipeDefaultSentinelDefersToTheSearch() {
        // FragPipe passes "r" when the user leaves the fragment-type selector alone, which is now
        // the case that configures itself from fragpipe.workflow.
        FragmentsExportMainClass.Options o = FragmentsExportMainClass.parseIonTypes("r");
        assertTrue(o.letters.isEmpty());
        assertFalse(o.neutralLoss);
        assertTrue(FragmentsExportMainClass.parseIonTypes("R").letters.isEmpty());
    }

    @Test
    void anEmptyArgumentAlsoDefersToTheSearch() {
        assertTrue(FragmentsExportMainClass.parseIonTypes("").letters.isEmpty());
        assertTrue(FragmentsExportMainClass.parseIonTypes("   ").letters.isEmpty());
        assertTrue(FragmentsExportMainClass.parseIonTypes(null).letters.isEmpty());
    }

    @Test
    void anExplicitListOverridesTheSearch() {
        assertEquals(List.of("b", "y"), FragmentsExportMainClass.parseIonTypes("b,y").letters);
        assertEquals(List.of("c", "z"), FragmentsExportMainClass.parseIonTypes("c_z").letters);
        assertEquals(List.of("b", "y"), FragmentsExportMainClass.parseIonTypes("by").letters);
        assertEquals(List.of("a", "b", "c"), FragmentsExportMainClass.parseIonTypes("a b c").letters);
    }

    @Test
    void neutralLossIsAnAdditiveFlag() {
        FragmentsExportMainClass.Options o = FragmentsExportMainClass.parseIonTypes("b_y_neu");
        assertEquals(List.of("b", "y"), o.letters);
        assertTrue(o.neutralLoss);
    }

    @Test
    void retiredGlycoKeywordsAreAcceptedAndIgnored() {
        // An older FragPipe still passes these. They must not fail, and must not be mistaken for
        // ion letters — glycan ions now come from the search's own labile parameters.
        FragmentsExportMainClass.Options o = FragmentsExportMainClass.parseIonTypes("b_y_ngly");
        assertEquals(List.of("b", "y"), o.letters);
        assertTrue(FragmentsExportMainClass.parseIonTypes("ngly_ogly_gly").letters.isEmpty());
    }

    @Test
    void unrecognizedTokensAreDroppedRatherThanGuessed() {
        assertEquals(List.of("b"), FragmentsExportMainClass.parseIonTypes("b,wat").letters);
    }
}
