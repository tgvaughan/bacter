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

    public int getIntervalCount() {
        return decendentLeaves.size();
    }

    /**
     * Computes the union between this ancestry and another, additionally
     * producing a SiteAncestry representing the those sites and samples
     * which experienced a coalescent event.
     *
     * @param other
     * @param coalescence
     * @param union
     */
    public void merge(SiteAncestry other, SiteAncestry coalescence, SiteAncestry union) {
        int i=0, j=0;

        while (i<getIntervalCount()) {

            boolean started = false;
            while (j<other.getIntervalCount()) {

                if (other.siteRanges.get(2*j) >= siteRanges.get(2*i+1))
                    break;

                if (other.siteRanges.get(2*j+1)<siteRanges.get(2*i)) {
                    union.siteRanges.add(other.siteRanges.get(2*j));
                    union.siteRanges.add(other.siteRanges.get(2*j+1));
                    union.decendentLeaves.add(other.decendentLeaves.get(j));

                    j += 1;
                    continue;
                }

                if (other.siteRanges.get(2*j)>siteRanges.get(2*i)) {
                    if (!started) {
                        union.siteRanges.add(siteRanges.get(2*i));
                        union.siteRanges.add(other.siteRanges.get(2*j));
                        union.decendentLeaves.add(decendentLeaves.get(i));
                        started = true;
                    }


                }


                j += 1;
            }

            i += 1;
        }

    }
}
