import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumMap;
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
    void customSeriesOffsetIsRelativeToTheBareResidueSum() {
        // MSFragger adds a custom offset to the residue sum at either terminus, in the same field
        // where its built-in z is +1.991841 and y is +18.010565. So zdot below is MSFragger's own
        // z-radical. Resolving a C-terminal offset against y instead moves every custom ion by
        // 18 Da, which annotates nothing and looks exactly like a series the search never used.
        List<CustomIon> custom = List.of(
                new CustomIon("zdot", false, 1.991841),
                new CustomIon("cdot", true, 0.02381));
        List<FragmentAnnotator.Series> resolved =
                FragmentAnnotator.resolveSeries(Arrays.asList("zdot", "cdot", "z"), custom);

        FragmentAnnotator.Series zdot = seriesNamed(resolved, "zdot");
        FragmentAnnotator.Series z = seriesNamed(resolved, "z");
        assertFalse(zdot.nterm);
        assertEquals(z.shift, zdot.shift, 1e-5,
                "zdot C 1.991841 is the z-radical, so its shift must match the standard z series");

        FragmentAnnotator.Series cdot = seriesNamed(resolved, "cdot");
        assertTrue(cdot.nterm);
        assertEquals(0.02381, cdot.shift, 1e-9,
                "b is the bare residue sum, so an N-terminal offset applies to it unchanged");
    }

    @Test
    void customCTerminalIonIsAnnotatedWhereMSFraggerScoredIt() {
        // "zstar C -16.01872" on a C-terminal K is scored by MSFragger at
        // 128.09496 - 16.01872 + proton = 113.0835, not at y1 - 16.01872 = 131.0941.
        FragmentAnnotator.configureTolerance(20, false);
        List<FragmentAnnotator.Series> backbone = FragmentAnnotator.resolveSeries(
                List.of("zstar"), List.of(new CustomIon("zstar", false, -16.01872)));
        double[] mzs = {113.08352, 131.09408};
        double[] ints = {100, 100};

        EnumMap<IonCategory, ArrayList<IonMatch>> result = FragmentAnnotator.annotate(
                "AK", new ArrayList<>(), 0.0, 1, mzs, ints, backbone, List.of(),
                SearchParams.parse(""), false, null);

        List<IonMatch> matched = result.get(IonCategory.BACKBONE);
        assertEquals(1, matched.size());
        assertEquals("zstar1", matched.get(0).label);
        assertEquals(113.08352, matched.get(0).theoMz, 1e-4);
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
