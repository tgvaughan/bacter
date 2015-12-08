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

    public SiteAncestry getUnionWith(SiteAncestry other) {
        SiteAncestry union = new SiteAncestry();

        int i=0, j=0;

        while (i<getIntervalCount()) {
            while (j<other.getIntervalCount()) {

                if (other.siteRanges.get(2*j) >= siteRanges.get(2*i+1))
                    break;

                if (other.siteRanges.get(2*j+1) > siteRanges.get(2*i)) {

                    if (other.siteRanges.get(2*j) < siteRanges.get(2*i)) {
                        if (other.decendentLeaves.get(j).equals(decendentLeaves.get(i))) {

                        } else {
                            
                        }
                    }

                } else {
                    union.siteRanges.add(other.siteRanges.get(2*j));
                    union.siteRanges.add(other.siteRanges.get(2*j+1));
                    union.decendentLeaves.add(other.decendentLeaves.get(j));
                }

                j += 1;
            }

            i += 1;
        }

        return union;
    }
}
