/**
 * A modification at a specific position in a peptide.
 * Position is 1-indexed (1 = first residue / N-terminus).
 * N-terminal mods: ptmName contains "N-term".
 * C-terminal mods: ptmName contains "C-term".
 */
public class ModificationMatch {
    private final String ptmName;
    private final int    position;
    private final double mass;

    public ModificationMatch(String ptmName, int position, double mass) {
        this.ptmName  = ptmName;
        this.position = position;
        this.mass     = mass;
    }

    public String getTheoreticPtm() { return ptmName; }
    public int getModificationSite() { return position; }
    public double getMass() { return mass; }

    public boolean isNTerminal() {
        return ptmName.contains("N-term") || ptmName.contains("n-term");
    }

    public boolean isCTerminal() {
        return ptmName.contains("C-term") || ptmName.contains("c-term");
    }
}
