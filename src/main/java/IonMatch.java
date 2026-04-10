/** A matched ion: theoretical m/z paired with an observed spectrum peak. */
public class IonMatch {
    /** Human-readable label, e.g. "b3", "y5++", "GLY_HN_HN", "SEQ_HN_HE" */
    public final String label;
    public final double peakMz;
    public final double peakIntensity;
    public final double theoMz;

    public IonMatch(String label, double peakMz, double peakIntensity, double theoMz) {
        this.label        = label;
        this.peakMz       = peakMz;
        this.peakIntensity = peakIntensity;
        this.theoMz       = theoMz;
    }

    public String getPeakAnnotation() {
        return label;
    }
}
