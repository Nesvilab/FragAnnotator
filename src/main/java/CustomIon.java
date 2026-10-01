/**
 * One user-defined fragment ion series from {@code msfragger.ion_series_definitions}, e.g.
 * {@code zdot C 1.991841;cdot N 0.02381}. Handled exactly like a standard series: generated
 * across the whole backbone at every charge, and reported in the backbone column under the name
 * the user gave it.
 *
 * <p>{@link #offset} is the number exactly as MSFragger uses it: added to the <b>bare residue
 * sum</b> of the fragment, at either terminus. That is the same field where MSFragger's built-in
 * b is 0, y is +18.010565 and z-radical is +1.991841, so {@code b* N -17.026548} is b-NH3 and
 * {@code zdot C 1.991841} is the z-radical (y-NH2). For an N-terminal series the residue sum is
 * the b neutral; for a C-terminal one it is y <b>without its water</b>. Getting the base wrong
 * moves every custom ion by 18 Da, which annotates nothing and looks exactly like a series the
 * search never used.
 */
public class CustomIon {

    public final String name;
    /** True for an N-terminal series (counts from the peptide N-term, like a/b/c). */
    public final boolean nterm;
    public final double offset;

    public CustomIon(String name, boolean nterm, double offset) {
        this.name = name;
        this.nterm = nterm;
        this.offset = offset;
    }

    /**
     * The shift from this series' base neutral mass: the b neutral for an N-terminal series, the y
     * neutral for a C-terminal one. MSFragger's offset is relative to the bare residue sum, which
     * is the b neutral already but is one water short of the y neutral, so a C-terminal offset
     * loses that water here.
     */
    public double shiftFromBase() {
        return nterm ? offset : offset - FragmentAnnotator.H2O;
    }

    @Override
    public String toString() {
        return name + (nterm ? " N " : " C ") + offset;
    }
}
