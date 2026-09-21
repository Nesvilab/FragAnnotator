/**
 * One user-defined fragment ion series from {@code msfragger.ion_series_definitions}, e.g.
 * {@code zOne C -16.01872;cdot N 0.02381}. Handled exactly like a standard series: generated
 * across the whole backbone at every charge, and reported in the backbone column under the name
 * the user gave it.
 *
 * <p>{@link #offset} is stated against the <b>ordinary ion of that terminus</b>, which is what
 * MSFragger's own parameter examples say: {@code b* N -17.026548} is b-NH3 and
 * {@code b0 N -18.010565} is b-H2O, so an N-terminal offset is relative to <b>b</b>;
 * {@code zOne C -16.01872} is the z-radical (y-NH2), so a C-terminal one is relative to <b>y</b>.
 * Getting the base wrong moves every custom ion by 18 Da, which annotates nothing and looks
 * exactly like a series the search never used.
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
     * neutral for a C-terminal one. Both bases are the ordinary ion of that terminus, and the y
     * neutral already carries its water, so the declared offset applies unchanged in both cases.
     */
    public double shiftFromBase() {
        return offset;
    }

    @Override
    public String toString() {
        return name + (nterm ? " N " : " C ") + offset;
    }
}
