package bacter;

import beast.evolution.tree.Node;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * Class of objects representing the site ancestry of a given lineage.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SiteAncestry {

    public List<Integer> siteRanges;
    public List<BitSet> decendentLeaves;

    public SiteAncestry() {
        siteRanges = new ArrayList<>();
        decendentLeaves = new ArrayList<>();
    }

    public SiteAncestry(Node node, Locus locus) {
        siteRanges = new ArrayList<>();
        siteRanges.add(0);
        siteRanges.add(locus.getSiteCount());

        decendentLeaves = new ArrayList<>();
        BitSet bitSet = new BitSet();
        bitSet.set(node.getNr());
        decendentLeaves.add(bitSet);
    }

    public SiteAncestry(String string) {

        siteRanges = new ArrayList<>();
        decendentLeaves = new ArrayList<>();

        string = string.replaceAll("\\s+","");

        if (string.isEmpty())
            return;

        string = string.substring(1, string.length()-1);
        String[] split1 = string.replaceAll("\\s","").split("}\\[");

        for (String aSplit1 : split1) {
            String[] split2 = aSplit1.split("]\\{");
            String[] rangeStr = split2[0].split(",");

            siteRanges.add(Integer.parseInt(rangeStr[0]));
            siteRanges.add(Integer.parseInt(rangeStr[1]));

            String[] bitStr = split2[1].split(",");
            BitSet theseDecendents = new BitSet();
            for (String aBitStr : bitStr) {
                theseDecendents.set(Integer.parseInt(aBitStr));
            }
            decendentLeaves.add(theseDecendents);
        }
    }

    public int getIntervalCount() {
        return decendentLeaves.size();
    }

    /**
     * Computes the union between this ancestry and another, additionally
     * producing a SiteAncestry representing the those sites and samples
     * which experienced a coalescent event.  Note that SAs are assumed
     * to come from contemporaneous points on extant lineages, meaning that
     * each site+sample combination must appear at most once.
     *
     * That is, the merge of A = [1,10]{0} and B = [5,20]{0,1} is not allowed
     * as this implies that sites [5,10] of sample 0 coexist in two
     * contemporaneoous lineages.
     *
     * @param other SA with which to merge
     * @param coalescence empty SA object in which to record coalescences
     * @param union empty SA object in which to record union
     */
    public void merge(SiteAncestry other, SiteAncestry coalescence, SiteAncestry union) {

        int i = 0, j = 0;

        while (i<getIntervalCount()) {

            while (j<other.getIntervalCount()
                    && other.siteRanges.get(2*j) < siteRanges.get(2*i)) {

                int xp = other.siteRanges.get(2*j);
                int yp = other.siteRanges.get(2*j + 1);
                yp = yp > siteRanges.get(2*i) ? siteRanges.get(2*i) : yp;

                union.siteRanges.add(xp);
                union.siteRanges.add(yp);
                union.decendentLeaves.add(other.decendentLeaves.get(j));

                if (other.siteRanges.get(2*j+1) <= siteRanges.get(2*i))
                    j += 1;
                else
                    break;
            }

            i += 1;
        }

    }

    /**
     * Constructs two new SiteAncestry objects describing the ancestry of
     * sites which fall respectively inside of and outside of the interval
     * [x,y].
     *
     * @param x left-hand boundary of interval
     * @param y right-hand boundary of interval
     * @param inside SA to hold ancestry of sites falling inside the interval
     * @param outside SA to hold ancestry of sites falling outside the interval
     */
    public void split(int x, int y, SiteAncestry inside, SiteAncestry outside) {

        int i=0;

        while (i<getIntervalCount() && siteRanges.get(2*i)<x) {
            int xp = siteRanges.get(2*i);
            int yp = siteRanges.get(2*i+1);
            yp = yp > x ? x : yp;

            outside.siteRanges.add(xp);
            outside.siteRanges.add(yp);
            outside.decendentLeaves.add(decendentLeaves.get(i));

            if (siteRanges.get(2*i+1) <= x)
                i += 1;
            else
                break;
        }

        while (i<getIntervalCount() && siteRanges.get(2*i)<y) {
            int xp = siteRanges.get(2*i);
            int yp = siteRanges.get(2*i+1);
            xp = xp < x ? x : xp;
            yp = yp > y ? y : yp;

            inside.siteRanges.add(xp);
            inside.siteRanges.add(yp);
            inside.decendentLeaves.add(decendentLeaves.get(i));

            if (siteRanges.get(2*i+1) <=y)
                i += 1;
            else
                break;
        }

        while (i<getIntervalCount()) {
            int xp = siteRanges.get(2*i);
            int yp = siteRanges.get(2*i+1);
            xp = xp < y ? y : xp;

            outside.siteRanges.add(xp);
            outside.siteRanges.add(yp);
            outside.decendentLeaves.add(decendentLeaves.get(i));

            i += 1;
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SiteAncestry that = (SiteAncestry) o;

        return siteRanges.equals(that.siteRanges)
                && decendentLeaves.equals(that.decendentLeaves);
    }

    @Override
    public int hashCode() {
        int result = siteRanges.hashCode();
        result = 31 * result + decendentLeaves.hashCode();
        return result;
    }

    @Override
    public String toString() {
        String res = "";

        for (int i=0; i<getIntervalCount(); i++) {
            if (i>0)
                res += " ";

            res += "[" + siteRanges.get(2*i) + "," + siteRanges.get(2*i+1) + "]"
                    + decendentLeaves.get(i).toString().replaceAll("\\s","");
        }

        return res;
    }

    /**
     * Main method for testing
     *
     * @param args unused
     */
    public static void main(String[] args) {

        SiteAncestry a = new SiteAncestry("[200,400]{0} [600,800]{1}");

        SiteAncestry inside = new SiteAncestry();
        SiteAncestry outside = new SiteAncestry();
        int x=0, y=100;
        a.split(x, y, inside, outside);

        System.out.println("a = " + a);
        System.out.println("inside a.split(" + x + "," + y + ") = " + inside);
        System.out.println("outside a.split(" + x + "," + y + ") = " + outside);
    }
}
