import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/** Pins ion-series resolution, remainder combination, and the labels remainders carry. */
class FragmentAnnotatorTest {

    private static FragmentAnnotator.Series seriesNamed(List<FragmentAnnotator.Series> all, String label) {
        for (FragmentAnnotator.Series s : all) if (s.label.equals(label)) return s;
        return null;
    }

    @Test
    void resolvesStandardSeriesByLetter() {
        List<FragmentAnnotator.Series> s =
                FragmentAnnotator.resolveSeries(Arrays.asList("b", "y", "c", "z"), List.of());
        assertEquals(4, s.size());
        assertTrue(seriesNamed(s, "b").nterm);
        assertFalse(seriesNamed(s, "y").nterm);
        assertEquals(0.0, seriesNamed(s, "b").shift, 1e-9);
        assertEquals(0.0, seriesNamed(s, "y").shift, 1e-9);
    }

    @Test
    void customSeriesShiftIsRelativeToItsTerminusOrdinaryIon() {
        // The reference search's zOne is the z-radical, i.e. y - NH2. Resolving it against the
        // wrong base moves every custom C-terminal ion by 18 Da, which annotates nothing and looks
        // exactly like a series the search never used.
        List<CustomIon> custom = List.of(
                new CustomIon("zOne", false, -16.01872),
                new CustomIon("cdot", true, 0.02381));
        List<FragmentAnnotator.Series> resolved =
                FragmentAnnotator.resolveSeries(Arrays.asList("zOne", "cdot", "z"), custom);

        FragmentAnnotator.Series zOne = seriesNamed(resolved, "zOne");
        FragmentAnnotator.Series z = seriesNamed(resolved, "z");
        assertFalse(zOne.nterm);
        assertEquals(z.shift, zOne.shift, 1e-5,
                "zOne is z-radical, so its shift must match the standard z series");

        FragmentAnnotator.Series cdot = seriesNamed(resolved, "cdot");
        assertTrue(cdot.nterm);
        assertEquals(0.02381, cdot.shift, 1e-9, "an N-terminal offset applies to b unchanged");
    }

    @Test
    void dropsSeriesTheResultNeverDeclared() {
        List<FragmentAnnotator.Series> s =
                FragmentAnnotator.resolveSeries(Arrays.asList("b", "nosuch"), List.of());
        assertEquals(1, s.size());
        assertEquals("b", s.get(0).label);
    }

    @Test
    void remainderTagIsSignedAndRounded() {
        assertEquals("+80", FragmentAnnotator.remainderTag(79.96633));
        assertEquals("+203", FragmentAnnotator.remainderTag(203.07937));
        assertEquals("-42", FragmentAnnotator.remainderTag(-42.0205));
        assertEquals("-18", FragmentAnnotator.remainderTag(-18.01056));
        assertEquals("+0", FragmentAnnotator.remainderTag(0.0));
    }

    @Test
    void oneModifiedSiteYieldsItsOwnRemainders() {
        double[] sums = FragmentAnnotator.uniqueSums(List.of(new double[]{0.0, 203.07937}));
        Arrays.sort(sums);
        assertArrayEquals(new double[]{0.0, 203.07937}, sums, 1e-6);
    }

    @Test
    void severalModificationsCombineAndCollapseOnMass() {
        // Two sites, each able to keep 0 or 114.03169. The cross-product has four members but only
        // three masses: 0+114 and 114+0 are one ion, and giving them two labels would put two
        // annotations on one peak.
        double[] sums = FragmentAnnotator.uniqueSums(List.of(
                new double[]{0.0, 114.03169},
                new double[]{0.0, 114.03169}));
        Arrays.sort(sums);
        assertArrayEquals(new double[]{0.0, 114.03169, 228.06338}, sums, 1e-6);
    }

    @Test
    void combinesRemaindersFromDifferentModifications() {
        double[] sums = FragmentAnnotator.uniqueSums(List.of(
                new double[]{0.0, 79.96633},
                new double[]{0.0, 203.07937}));
        Arrays.sort(sums);
        assertArrayEquals(new double[]{0.0, 79.96633, 203.07937, 283.0457}, sums, 1e-5);
    }

    @Test
    void noModificationsYieldsTheEmptySum() {
        assertArrayEquals(new double[]{0.0}, FragmentAnnotator.uniqueSums(new ArrayList<>()), 1e-9);
    }
}
