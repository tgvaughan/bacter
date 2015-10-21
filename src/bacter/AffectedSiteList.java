package bacter;

import beast.evolution.tree.Node;
import com.google.common.collect.Range;

import java.util.*;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class AffectedSiteList {

    ConversionGraph acg;
    public Map<Conversion, List<Integer>> affectedSites;
    public Map<Conversion, Integer> affectedSiteCount;
    public Map<Conversion, Double> affectedSiteFraction;

    ACGEventList acgEventList;

    public AffectedSiteList(ConversionGraph acg) {
        this.acg = acg;

        acgEventList = new ACGEventList(acg);
        affectedSites = new HashMap<>();
        affectedSiteCount = new HashMap<>();
        affectedSiteFraction = new HashMap<>();

        computeAffectedSites();
    }

    public Map<Locus, List<Integer>> getLeafAncestralSites() {
        Map<Locus, List<Integer>> res = new HashMap<>();

        for (Locus locus : acg.getLoci()) {
            List<Integer> siteRange = new ArrayList<>();
            siteRange.add(0);
            siteRange.add(locus.getSiteCount()-1);
            res.put(locus, siteRange);
        }

        return res;
    }

    public void computeAffectedSites() {

        Map<Node, Map<Locus, List<Integer>>> activeCFNodes = new HashMap<>();
        Map<Locus, Set<Conversion>> activeConversions = new HashMap<>();
        for (Locus locus : acg.getLoci())
            activeConversions.put(locus, new HashSet<>());

        Map<Locus, List<Integer>> ancestralSitesCF;
        List<Integer> ancestralSitesLocusCF;

        int leavesSeen = 0;
        boolean mrcaReached = false;
        for (ACGEventList.Event event : acgEventList.getACGEvents()) {

            if (mrcaReached) {
                if (event.type == ACGEventList.EventType.CONV_DEPART) {
                    affectedSites.put(event.conversion, new ArrayList<>());
                    affectedSiteCount.put(event.conversion, 0);
                    affectedSiteFraction.put(event.conversion, 0.0);
                }

                continue;
            }

            switch(event.type) {
                case CF_LEAF:
                    activeCFNodes.put(event.node, getLeafAncestralSites());
                    leavesSeen += 1;
                    break;

                case CF_COALESCENCE:
                    Node node1 = event.node.getLeft();
                    Node node2 = event.node.getRight();

                    ancestralSitesCF = new HashMap<>();
                    for (Locus locus : acg.getLoci()) {
                        ancestralSitesCF.put(locus,
                                getUnion(activeCFNodes.get(node1).get(locus),
                                        activeCFNodes.get(node2).get(locus)));
                    }

                    activeCFNodes.remove(node1);
                    activeCFNodes.remove(node2);
                    activeCFNodes.put(event.node, ancestralSitesCF);

                    break;

                case CONV_DEPART:
                    List<Integer> inside = new ArrayList<>();
                    List<Integer> outside = new ArrayList<>();
                    partitionRanges(activeCFNodes.get(event.node).get(event.conversion.getLocus()),
                            event.conversion.getStartSite(),
                            event.conversion.getEndSite()+1,
                            inside, outside);

                    affectedSites.put(event.conversion, inside);
                    affectedSiteCount.put(event.conversion, getTotalSites(inside));
                    affectedSiteFraction.put(event.conversion,
                            getTotalSites(inside) / (double) event.conversion.getSiteCount());
                    activeCFNodes.get(event.node).put(
                            event.conversion.getLocus(), outside);
                    activeConversions.get(event.conversion.locus).add(event.conversion);

                    if (leavesSeen == acg.getLeafNodeCount() && haveReachedAllMRCAs(activeCFNodes, activeConversions))
                        mrcaReached = true;

                    break;

                case CONV_ARRIVE:
                    activeCFNodes.get(event.node).put(event.conversion.locus,
                            getUnion(affectedSites.get(event.conversion),
                                    activeCFNodes.get(event.node).get(event.conversion.locus)));
                    activeConversions.get(event.conversion.locus).remove(event.conversion);

                    if (leavesSeen == acg.getLeafNodeCount() && haveReachedAllMRCAs(activeCFNodes, activeConversions))
                        mrcaReached = true;
                    break;
            }

        }

    }

    /**
     * Obtain the total number of sites included in a list of site ranges.
     *
     * @param as affected sites
     * @return total number of sites included
     */
    static int getTotalSites(List<Integer> as) {
        int res = 0;

        for (int i=0; i<as.size(); i+=2) {
            res += as.get(i+1)-as.get(i);
        }

        return res;
    }

    boolean haveReachedAllMRCAs(Map<Node, Map<Locus, List<Integer>>> activeCFNodes,
                                Map<Locus, Set<Conversion>> activeConversions) {

        for (Locus locus : acg.getLoci()) {
            List<Integer> startSites = new ArrayList<>();
            List<Integer> endSites = new ArrayList<>();
            for (Node node : activeCFNodes.keySet()) {
                for (int i=0; i<activeCFNodes.get(node).get(locus).size(); i+=2) {
                    startSites.add(activeCFNodes.get(node).get(locus).get(i));
                    endSites.add(activeCFNodes.get(node).get(locus).get(i + 1));
                }
            }

            for (Conversion conv : activeConversions.get(locus)) {
                for (int i=0; i<affectedSites.get(conv).size(); i+=2) {
                    startSites.add(affectedSites.get(conv).get(i));
                    endSites.add(affectedSites.get(conv).get(i + 1));
                }
            }

            Collections.sort(startSites);
            Collections.sort(endSites);

            for (int i=0; i<startSites.size()-1; i++) {
                if (startSites.get(i+1)<endSites.get(i))
                    return false;
            }

            return true;
        }

        return true;
    }

    /**
     * Test whether two range lists are disjoint.
     *
     * @param as1 list 1
     * @param as2 list 2
     * @return true iff all ranges are disjoint
     */
    static boolean rangesAreDisjoint(List<Integer> as1, List<Integer> as2) {

        List<Integer> startSites = new ArrayList<>();
        List<Integer> endSites = new ArrayList<>();

        for (int i=0; i<as1.size(); i+=2) {
            startSites.add(as1.get(i));
            endSites.add(as1.get(i + 1));
        }

        for (int i=0; i<as2.size(); i+=2) {
            startSites.add(as2.get(i));
            endSites.add(as2.get(i+1));
        }

        Collections.sort(startSites);
        Collections.sort(endSites);

        for (int i=0; i<startSites.size()-1; i++) {
            if (startSites.get(i+1)<endSites.get(i))
                return false;
        }

        return true;
    }

    static List<Integer> getUnion(List<Integer> as1, List<Integer> as2) {
        List<Integer> union = new ArrayList<>();
        int nextx, nexty;

        int i1 = 0, i2 = 0, ui = -2;
        int N1 = as1.size(), N2 = as2.size();
        while (i1 < N1 || i2 < N2) {
            if (i1<N1) {
                if (i2==N2 || as1.get(i1) < as2.get(i2)) {
                    nextx = as1.get(i1);
                    nexty = as1.get(i1+1);
                    i1 += 2;
                } else {
                    nextx = as2.get(i2);
                    nexty = as2.get(i2+1);
                    i2 += 2;
                }
            } else {
                nextx = as2.get(i2);
                nexty = as2.get(i2+1);
                i2 += 2;
            }

            if (ui<0 || union.get(ui+1)<nextx) {
                union.add(nextx);
                union.add(nexty);
                ui += 2;
            } else {
                if (union.get(ui+1)<nexty)
                    union.set(ui+1, nexty);
            }
        }

        return union;
    }

    public static void partitionRanges(List<Integer> as, int x, int y,
                                       List<Integer> inside, List<Integer> outside) {

        int i=0;
        while (i<as.size() && as.get(i) < x)
            outside.add(as.get(i++));

        if (i%2==1) {
            outside.add(x);
            if (x<as.get(i))
                inside.add(x);
            else
                i += 1;
        }

        while (i<as.size() && as.get(i)<y)
            inside.add(as.get(i++));

        if (i%2==1) {
            inside.add(y);
            if (y<as.get(i))
                outside.add(y);
            else
                i += 1;
        }

        while (i<as.size())
            outside.add(as.get(i++));

    }

    static String rangesToString(List<Integer> bounds) {
        String res = "";

        for (int i=0; i<bounds.size(); i+=2)
            res += " [" + bounds.get(i) + "," + bounds.get(i+1) + "]";

        return "{" + res + " }";
    }

    public static void main(String[] args) {
        List<Integer> as1 = new ArrayList<>();
        as1.add(3); as1.add(5);
        as1.add(9); as1.add(17);
        System.out.println("as1 = " + rangesToString(as1));

        List<Integer> as2 = new ArrayList<>();
        as2.add(1); as2.add(2);
        as2.add(5); as2.add(10);
        System.out.println("as2 = " + rangesToString(as2));
        List<Integer> as3 = getUnion(as1, as2);

        System.out.println("as3 = union(as1,as2) = " + rangesToString(as3));

        List<Integer> intersect = new ArrayList<>();
        List<Integer> difference = new ArrayList<>();
        int x = 4, y = 24;
        partitionRanges(as3, x, y, intersect, difference);
        System.out.format("as3 - [%d,%d] = %s\n", x, y, rangesToString(difference));
        System.out.format("as3 ^ [%d,%d] = %s\n", x, y, rangesToString(intersect));

        List<Integer> as4 = new ArrayList<>();
        as4.add(10); as4.add(20);
        as4.add(30); as4.add(40);
        List<Integer> as5 = new ArrayList<>();
        as5.add(1); as5.add(11);
        as5.add(20); as5.add(30);
        System.out.println(rangesAreDisjoint(as4, as5));

    }

}
