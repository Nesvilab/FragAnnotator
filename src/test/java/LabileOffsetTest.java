import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/** Pins one mass-offset entry: its sites, its declared ion lists, and the implicit zero. */
class LabileOffsetTest {

    @Test
    void parsesAPlainResidueList() {
        LabileOffset o = LabileOffset.parseEntry("541.06110(aa=SKTYHDE_d=136.06232_p=0.0,114.03169_f=0.0)");
        assertNotNull(o);
        assertEquals(541.0611, o.mass, 1e-6);
        assertTrue(o.allows('S'));
        assertTrue(o.allows('T'));
        assertTrue(o.allows('e'), "site matching ignores case");
        assertFalse(o.allows('R'));
    }

    @Test
    void takesOnlyTheModifiedResidueOfASequon() {
        // {N[^P][ST]} is the N-glycan sequon. Only its first residue carries the modification;
        // the rest was the search's own constraint, already satisfied by the assignment.
        LabileOffset o = LabileOffset.parseEntry("203.07937(aa={N[^P][ST]}_d=204.08665_f=203.07937)");
        assertNotNull(o);
        assertTrue(o.allows('N'));
        assertFalse(o.allows('S'), "S appears in the motif but is not the modified residue");
        assertFalse(o.allows('T'));
        assertFalse(o.allows('P'));
    }

    @Test
    void parenthesesNameTheModifiedResidue() {
        // A sequon may mark its modified residue explicitly, and then it is not the first one.
        LabileOffset o = LabileOffset.parseEntry("100.0(aa={[ST]P(N)}_d=50.0)");
        assertNotNull(o);
        assertTrue(o.allows('N'));
        assertFalse(o.allows('S'));
        assertFalse(o.allows('P'));
    }

    @Test
    void closesTheEntryOnItsLastParenthesis() {
        // A parenthesised sequon must not truncate the body: everything after the first ')'
        // would otherwise be lost, taking the declared ion lists with it.
        LabileOffset o = LabileOffset.parseEntry("100.0(aa={(N)[^P][ST]}_d=204.08665,186.07609_p=203.07937)");
        assertNotNull(o);
        assertTrue(o.allows('N'));
        assertEquals(2, o.diagnostic.length);
        assertArrayEquals(new double[]{0.0, 203.07937}, o.peptide, 1e-6);
    }

    @Test
    void anEmptySiteListMeansAnyResidue() {
        LabileOffset o = LabileOffset.parseEntry("100.0(aa=_d=50.0)");
        assertNotNull(o);
        assertTrue(o.allows('A'));
        assertTrue(o.allows('W'));
    }

    @Test
    void aWildcardMeansAnyResidue() {
        LabileOffset o = LabileOffset.parseEntry("100.0(aa=*_d=50.0)");
        assertNotNull(o);
        assertTrue(o.allows('A'));
    }

    @Test
    void aNegatedClassExcludesItsResidues() {
        LabileOffset o = LabileOffset.parseEntry("100.0(aa={([^P])XX}_d=50.0)");
        assertNotNull(o);
        assertFalse(o.allows('P'));
        assertTrue(o.allows('N'));
    }

    @Test
    void addsZeroToBothRemainderLists() {
        // Complete loss is available to any labile modification, and the declared lists routinely
        // omit it: a glyco offset declares only _f=203.07937, so without this the backbone ions of
        // a glycopeptide with the glycan fully stripped would never be annotated.
        LabileOffset o = LabileOffset.parseEntry("203.07937(aa={N[^P][ST]}_d=204.08665_p=203.07937_f=203.07937)");
        assertNotNull(o);
        assertArrayEquals(new double[]{0.0, 203.07937}, o.fragment, 1e-6);
        assertArrayEquals(new double[]{0.0, 203.07937}, o.peptide, 1e-6);
        // Diagnostic ions are observed m/z values, not remainders, so no zero is added.
        assertArrayEquals(new double[]{204.08665}, o.diagnostic, 1e-6);
    }

    @Test
    void doesNotDuplicateADeclaredZero() {
        LabileOffset o = LabileOffset.parseEntry("541.0611(aa=R_f=0.00000,-42.02050)");
        assertNotNull(o);
        assertArrayEquals(new double[]{-42.0205, 0.0}, o.fragment, 1e-6);
    }

    @Test
    void anOffsetDeclaringNoIonsIsNotLabile() {
        // The mandatory 0.0000(aa=) entry, and any offset with no labile fragmentation. The
        // implicit zero must not make an inert offset look labile.
        assertNull(LabileOffset.parseEntry("0.0000(aa=)"));
        assertNull(LabileOffset.parseEntry("541.0611(aa=R)"));
        assertNull(LabileOffset.parseEntry(""));
        assertNull(LabileOffset.parseEntry("   "));
        assertNull(LabileOffset.parseEntry("notanumber(aa=R_d=50.0)"));
    }

    @Test
    void readsAllThreeIonLists() {
        LabileOffset o = LabileOffset.parseEntry(
                "541.06110(aa=SKTYHDE_d=136.06232,250.09401_p=0.00000,114.03169_f=0.00000)");
        assertNotNull(o);
        assertArrayEquals(new double[]{136.06232, 250.09401}, o.diagnostic, 1e-6);
        assertArrayEquals(new double[]{0.0, 114.03169}, o.peptide, 1e-6);
        assertArrayEquals(new double[]{0.0}, o.fragment, 1e-6);
    }
}
