package bacter;

import beast.base.evolution.tree.Node;

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
    public List<BitSet> descendantLeaves;

    public SiteAncestry() {
        siteRanges = new ArrayList<>();
        descendantLeaves = new ArrayList<>();
    }

    public SiteAncestry(Node node, Locus locus) {
        siteRanges = new ArrayList<>();
        siteRanges.add(0);
        siteRanges.add(locus.getSiteCount());

        descendantLeaves = new ArrayList<>();
        BitSet bitSet = new BitSet();
        bitSet.set(node.getNr());
        descendantLeaves.add(bitSet);
    }

    public SiteAncestry(String string) {

        siteRanges = new ArrayList<>();
        descendantLeaves = new ArrayList<>();

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
            descendantLeaves.add(theseDecendents);
        }
    }

    public int getIntervalCount() {
        return descendantLeaves.size();
    }

    public void addInterval(int x, int y, BitSet dl) {
        if (getIntervalCount()>0
                && siteRanges.get(siteRanges.size()-1)==x
                && descendantLeaves.get(descendantLeaves.size()-1).equals(dl))
            siteRanges.set(siteRanges.size()-1, y);
        else {
            siteRanges.add(x);
            siteRanges.add(y);
            descendantLeaves.add(dl);
        }
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
    public void merge(SiteAncestry other, Coalescence coalescence, SiteAncestry union) {

        int i = 0, j = 0;

        int last=0;

        while (i<getIntervalCount()) {

            int x = siteRanges.get(2*i);
            int y = siteRanges.get(2*i + 1);
            BitSet dl = descendantLeaves.get(i);

            while (j<other.getIntervalCount()
                    && other.siteRanges.get(2*j) < x) {

                int xp = other.siteRanges.get(2*j);
                int yp = other.siteRanges.get(2*j + 1);
                xp = xp < last ? last : xp;
                yp = yp > x ? x : yp;

                union.addInterval(xp, yp, other.descendantLeaves.get(j));

                if (other.siteRanges.get(2*j+1) <= x)
                    j += 1;
                else
                    break;
            }

            last = x;

            while (j<other.getIntervalCount()
                    && other.siteRanges.get(2*j) < y) {

                int xp = other.siteRanges.get(2*j);
                int yp = other.siteRanges.get(2*j + 1);
                xp = xp < x ? x : xp;
                yp = yp > y ? y : yp;

                if (last < xp)
                    union.addInterval(last, xp, dl);

                last = yp;

                BitSet mergedDescendants = (BitSet)dl.clone();
                mergedDescendants.or(other.descendantLeaves.get(j));
                union.addInterval(xp, yp, mergedDescendants);

                coalescence.addInterval(xp, yp, dl, other.descendantLeaves.get(j));

                if (other.siteRanges.get(2*j+1) <= y)
                    j += 1;
                else
                    break;
            }

            if (last < y)
                union.addInterval(last, y, dl);

            last = y;

            i += 1;
        }

        while (j < other.getIntervalCount()) {

            int xp = other.siteRanges.get(2*j);
            if (!siteRanges.isEmpty() && xp < siteRanges.get(siteRanges.size()-1))
                xp = siteRanges.get(siteRanges.size()-1);
            int yp = other.siteRanges.get(2*j + 1);

            union.addInterval(xp, yp, other.descendantLeaves.get(j));

            j += 1;
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
            outside.descendantLeaves.add(descendantLeaves.get(i));

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
            inside.descendantLeaves.add(descendantLeaves.get(i));

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
            outside.descendantLeaves.add(descendantLeaves.get(i));

            i += 1;
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SiteAncestry that = (SiteAncestry) o;

        return siteRanges.equals(that.siteRanges)
                && descendantLeaves.equals(that.descendantLeaves);
    }

    @Override
    public int hashCode() {
        int result = siteRanges.hashCode();
        result = 31 * result + descendantLeaves.hashCode();
        return result;
    }

    @Override
    public String toString() {
        String res = "";

        for (int i=0; i<getIntervalCount(); i++) {
            if (i>0)
                res += " ";

            res += "[" + siteRanges.get(2*i) + "," + siteRanges.get(2*i+1) + "]"
                    + descendantLeaves.get(i).toString().replaceAll("\\s","");
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

        SiteAncestry b = new SiteAncestry("[0,100]{0} [200,300]{1}");
        SiteAncestry c = new SiteAncestry("[100,200]{0} [350,400]{1}");

        SiteAncestry union = new SiteAncestry();
        Coalescence coalescence = new Coalescence();
        b.merge(c, coalescence, union);

        System.out.println("\nb = " + b);
        System.out.println("c = " + c);
        System.out.println("b.merge(c) union = " + union);
        System.out.println("b.merge(c) coalescence = " + coalescence);

        SiteAncestry d = new SiteAncestry("[120,150]{0} [200,300]{1}");
        SiteAncestry e = new SiteAncestry("[100,200]{1}");

        union = new SiteAncestry();
        coalescence = new Coalescence();
        d.merge(e, coalescence, union);

        System.out.println("\nd = " + d);
        System.out.println("e = " + e);
        System.out.println("d.merge(e) union = " + union);
        System.out.println("d.merge(e) coalescence = " + coalescence);

    }
}
