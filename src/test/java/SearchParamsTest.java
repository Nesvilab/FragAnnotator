import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/** Pins the fragpipe.workflow parameter format the annotator configures itself from. */
class SearchParamsTest {

    /**
     * The detailed-offset parameter as FragPipe writes it, with Java properties escaping. This is
     * the ADPr reference search verbatim: the same mass declared twice with different sites and
     * different remainder ions.
     */
    private static final String ADPR_DETAILED =
            "msfragger.mass_offsets_detailed=0.00000(aa\\=);"
            + "541.06110(aa\\=SKTYHDE_d\\=136.06232,250.09401_p\\=0.00000,114.03169_f\\=0.00000);"
            + "541.06110(aa\\=R_d\\=584.09021_p\\=0.00000,114.03169_f\\=-42.02050)\n";

    @Test
    void unescapesJavaPropertiesEscaping() {
        // Without this every number after a "\=" keeps its backslash and parses as nothing.
        assertEquals("aa=SKTYHDE_d=136.06", SearchParams.unescape("aa\\=SKTYHDE_d\\=136.06"));
        assertEquals("a:b", SearchParams.unescape("a\\:b"));
        assertEquals("a\\b", SearchParams.unescape("a\\\\b"));
        assertEquals("a\tb", SearchParams.unescape("a\\tb"));
        assertEquals("a", SearchParams.unescape("a\\")); // a trailing backslash escapes nothing
    }

    @Test
    void tokenizesEverySeparatorFragPipeHasUsed() {
        assertEquals(List.of("b", "y"), SearchParams.tokens("b,y"));
        assertEquals(List.of("0", "114.03169", "193.99802"), SearchParams.tokens("0/114.03169/193.99802"));
        assertEquals(List.of("0", "203.07937"), SearchParams.tokens("0 203.07937"));
        assertEquals(List.of("zOne", "C", "-16.01872"), SearchParams.tokens("zOne;C;-16.01872"));
        assertTrue(SearchParams.tokens("").isEmpty());
    }

    @Test
    void readsIonSeriesPreservingCase() {
        // Uppercase Y is MSFragger's peptide-remainder series, not the backbone y. Lowercasing
        // here would silently merge the two.
        SearchParams s = SearchParams.parse("msfragger.fragment_ion_series=b,y,Y\n");
        assertEquals(List.of("b", "y", "Y"), s.ionSeries);
    }

    @Test
    void readsToleranceAndUnits() {
        SearchParams ppm = SearchParams.parse(
                "msfragger.fragment_mass_tolerance=7\nmsfragger.fragment_mass_units=1\n");
        assertEquals(7.0, ppm.fragTol, 1e-9);
        assertFalse(ppm.fragTolDa);

        SearchParams da = SearchParams.parse(
                "msfragger.fragment_mass_tolerance=0.02\nmsfragger.fragment_mass_units=0\n");
        assertEquals(0.02, da.fragTol, 1e-9);
        assertTrue(da.fragTolDa);
    }

    @Test
    void missingToleranceFallsBackToTwentyPpm() {
        // 0 is a tolerance, so absence cannot be signalled by a value.
        SearchParams s = SearchParams.parse("msfragger.fragment_ion_series=b,y\n");
        assertEquals(SearchParams.DEFAULT_TOL_PPM, s.fragTol, 1e-9);
        assertFalse(s.fragTolDa);
    }

    @Test
    void parsesCustomIonDefinitions() {
        SearchParams s = SearchParams.parse(
                "msfragger.ion_series_definitions=zOne C -16.01872;zTwo C -15.0109;cdot N 0.02381\n");
        assertEquals(3, s.customIons.size());
        assertEquals("zOne", s.customIons.get(0).name);
        assertFalse(s.customIons.get(0).nterm);
        assertEquals(-16.01872, s.customIons.get(0).offset, 1e-9);
        assertTrue(s.customIons.get(2).nterm);
    }

    @Test
    void dropsCustomDefinitionsWithAnUnknownTerminus() {
        SearchParams s = SearchParams.parse("msfragger.ion_series_definitions=bad X 1.0;ok N 2.0\n");
        assertEquals(1, s.customIons.size());
        assertEquals("ok", s.customIons.get(0).name);
    }

    @Test
    void labileModeIsOffWhenTheParameterSaysSo() {
        // The custom-ion reference folder is exactly this case: a populated detailed list left over
        // from a search that ran with labile mode off, whose offsets were never searched.
        SearchParams s = SearchParams.parse(
                "msfragger.labile_search_mode=off\nmsfragger.use_detailed_offsets=false\n" + ADPR_DETAILED);
        assertFalse(s.labile);
        assertTrue(s.offsets.isEmpty());
    }

    @Test
    void readsDetailedOffsetsForALabileSearch() {
        SearchParams s = SearchParams.parse(
                "msfragger.labile_search_mode=labile\nmsfragger.use_detailed_offsets=true\n" + ADPR_DETAILED);
        assertTrue(s.labile);
        // The mandatory 0.0000(aa=) entry declares no ions and is not an offset.
        assertEquals(2, s.offsets.size());
        assertTrue(s.hasDiagnostic());
        assertTrue(s.hasPepRemainder());
    }

    @Test
    void picksTheOffsetMatchingBothMassAndSite() {
        SearchParams s = SearchParams.parse(
                "msfragger.labile_search_mode=labile\nmsfragger.use_detailed_offsets=true\n" + ADPR_DETAILED);

        // Same mass, different residue, different declared ions. Matching on mass alone would
        // annotate the wrong remainder for the residue at hand.
        LabileOffset onS = s.offsetFor(541.0611, 'S');
        LabileOffset onR = s.offsetFor(541.0611, 'R');
        assertNotNull(onS);
        assertNotNull(onR);
        assertNotSame(onS, onR);
        assertEquals(584.09021, onR.diagnostic[0], 1e-6);
        assertArrayEquals(new double[]{-42.02050, 0.0}, onR.fragment, 1e-6);

        assertNull(s.offsetFor(541.0611, 'W'), "W is a site of neither offset");
        assertNull(s.offsetFor(15.9949, 'M'), "an ordinary variable mod matches no offset");
    }

    @Test
    void anUnlocalizedMassSkipsTheSiteCheck() {
        SearchParams s = SearchParams.parse(
                "msfragger.labile_search_mode=labile\nmsfragger.use_detailed_offsets=true\n" + ADPR_DETAILED);
        // Residue 0 means "the search could not place it", so there is no residue to test.
        assertNotNull(s.offsetFor(541.0611, (char) 0));
    }

    @Test
    void acceptsBothSpellingsOfTheDetailedOffsetKey() {
        String body = "0.0000(aa\\=);100.0(aa\\=S_d\\=50.0)\n";
        SearchParams a = SearchParams.parse("msfragger.labile_search_mode=labile\n"
                + "msfragger.use_detailed_offsets=true\nmsfragger.mass_offsets_detailed=" + body);
        SearchParams b = SearchParams.parse("msfragger.labile_search_mode=labile\n"
                + "msfragger.use_detailed_offsets=true\nmsfragger.detailed_mass_offsets=" + body);
        assertEquals(1, a.offsets.size());
        assertEquals(1, b.offsets.size());
    }

    @Test
    void fallsBackToGlobalListsWhenDetailedOffsetsAreOff() {
        // A labile search that did not use the detailed list would otherwise produce three empty
        // columns — present, so the file looks annotated.
        SearchParams s = SearchParams.parse(
                "msfragger.labile_search_mode=labile\n"
                + "msfragger.use_detailed_offsets=false\n"
                + "msfragger.mass_offsets=0 541.06111\n"
                + "msfragger.Y_type_masses=0 114.03169\n"
                + "msfragger.diagnostic_fragments=136.06232/250.09401\n"
                + "msfragger.remainder_fragment_masses=-42.0205\n");
        assertTrue(s.labile);
        assertEquals(1, s.offsets.size(), "the 0 entry in mass_offsets is not an offset");
        LabileOffset o = s.offsets.get(0);
        assertEquals(541.06111, o.mass, 1e-6);
        assertEquals(2, o.diagnostic.length);
        assertTrue(o.allows('R'), "a global offset has no site restriction");
        assertTrue(o.allows('S'));
    }

    @Test
    void aGlycoSearchKeepsItsGlobalListsEvenWithNoPerMassOffset() {
        // A glyco search keeps every glycan mass in mass_offsets_detailed and writes
        // mass_offsets=0, so with detailed offsets switched off there is no per-mass entry to
        // build. Without the global offset the glycopeptide would get no labile ions at all -
        // three empty columns on a file that looks annotated.
        SearchParams s = SearchParams.parse(
                "msfragger.labile_search_mode=nglycan\n"
                + "msfragger.use_detailed_offsets=false\n"
                + "msfragger.mass_offsets=0\n"
                + "msfragger.Y_type_masses=0 203.07937 406.15874\n"
                + "msfragger.diagnostic_fragments=204.086646 186.076086\n"
                + "msfragger.remainder_fragment_masses=203.07937\n");
        assertTrue(s.labile);
        assertTrue(s.offsets.isEmpty(), "mass_offsets holds no mass to build an offset from");
        assertNotNull(s.globalOffset, "the global ion lists must survive for the glycan to use");
        assertEquals(2, s.globalOffset.diagnostic.length);
        assertArrayEquals(new double[]{0.0, 203.07937}, s.globalOffset.fragment, 1e-6);
        assertTrue(s.globalOffset.allows('N'), "a global offset has no site restriction");
    }

    @Test
    void thereIsNoGlobalOffsetWhenTheSearchDeclaresNoLabileIons() {
        SearchParams s = SearchParams.parse(
                "msfragger.labile_search_mode=labile\nmsfragger.use_detailed_offsets=false\n"
                + "msfragger.mass_offsets=0\n");
        assertNull(s.globalOffset);
    }

    @Test
    void ignoresCommentsAndBlankLines() {
        SearchParams s = SearchParams.parse(
                "# MSFragger version 4.5\n\n! another comment\nmsfragger.fragment_ion_series=c,z\n");
        assertEquals(List.of("c", "z"), s.ionSeries);
    }

    @Test
    void nglycanModeIsRecognisedForWarningPurposesOnly() {
        SearchParams s = SearchParams.parse("msfragger.labile_search_mode=nglycan\n");
        assertTrue(s.labile);
        assertTrue(s.nglycanMode);
        // O-glyco runs in plain "labile" mode, which is why nothing branches on this flag.
        assertFalse(SearchParams.parse("msfragger.labile_search_mode=labile\n").nglycanMode);
    }
}
