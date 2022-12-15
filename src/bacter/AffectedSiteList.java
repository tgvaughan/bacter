package bacter;

import bacter.util.IntRanges;
import beast.base.evolution.tree.Node;

import java.util.*;

/**
 * Class of objects which contain the list of sites whose ancestry involves
 * each conversion. Used to determine which conversions to actually consider
 * when evaluating the ARG likelihood.
 *
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

        Map<Node, Map<Locus, List<Integer>>> activeCFNodes = new HashMap<>();
        Map<Locus, Set<Conversion>> activeConversions = new HashMap<>();
        for (Locus locus : acg.getConvertibleLoci())
            activeConversions.put(locus, new HashSet<>());

        Map<Locus, List<Integer>> ancestralSitesCF;

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

            switch (event.type) {
                case CF_LEAF:
                    activeCFNodes.put(event.node, getLeafAncestralSites());
                    leavesSeen += 1;
                    break;

                case CF_COALESCENCE:
                    Node node1 = event.node.getLeft();
                    Node node2 = event.node.getRight();

                    ancestralSitesCF = new HashMap<>();
                    for (Locus locus : acg.getConvertibleLoci()) {
                        ancestralSitesCF.put(locus,
                                IntRanges.getUnion(activeCFNodes.get(node1).get(locus),
                                        activeCFNodes.get(node2).get(locus)));
                    }

                    activeCFNodes.remove(node1);
                    activeCFNodes.remove(node2);
                    activeCFNodes.put(event.node, ancestralSitesCF);

                    if (leavesSeen == acg.getLeafNodeCount() && haveReachedAllMRCAs(activeCFNodes, activeConversions))
                        mrcaReached = true;

                    break;

                case CONV_DEPART:
                    List<Integer> inside = new ArrayList<>();
                    List<Integer> outside = new ArrayList<>();
                    IntRanges.partitionRanges(activeCFNodes.get(event.node).get(event.conversion.getLocus()),
                            event.conversion.getStartSite(),
                            event.conversion.getEndSite() + 1,
                            inside, outside);

                    affectedSites.put(event.conversion, inside);
                    affectedSiteCount.put(event.conversion, IntRanges.getTotalSites(inside));
                    affectedSiteFraction.put(event.conversion,
                            IntRanges.getTotalSites(inside) / (double) event.conversion.getSiteCount());
                    activeCFNodes.get(event.node).put(
                            event.conversion.getLocus(), outside);
                    activeConversions.get(event.conversion.locus).add(event.conversion);

                    break;

                case CONV_ARRIVE:
                    activeCFNodes.get(event.node).put(event.conversion.getLocus(),
                            IntRanges.getUnion(affectedSites.get(event.conversion),
                                    activeCFNodes.get(event.node).get(event.conversion.getLocus())));
                    activeConversions.get(event.conversion.getLocus()).remove(event.conversion);

                    if (leavesSeen == acg.getLeafNodeCount() && haveReachedAllMRCAs(activeCFNodes, activeConversions))
                        mrcaReached = true;
                    break;
            }

        }
    }

    /**
     * Assembles complete site list for association with a leaf node.
     *
     * @return list of sites
     */
    protected Map<Locus, List<Integer>> getLeafAncestralSites() {
        Map<Locus, List<Integer>> res = new HashMap<>();

        for (Locus locus : acg.getConvertibleLoci()) {
            List<Integer> siteRange = new ArrayList<>();
            siteRange.add(0);
            siteRange.add(locus.getSiteCount() - 1);
            res.put(locus, siteRange);
        }

        return res;
    }

    /**
     * Test to see whether MRCA of every site has been found.
     * This is actually pretty expensive.  There's got to be a better way...
     *
     * @param activeCFNodes set of active CF nodes and the sites they represent
     * @param activeConversions set of active conversions and the sites they represent
     * @return true if all sites have found an MRCA, false otherwise
     */
    protected boolean haveReachedAllMRCAs(Map<Node, Map<Locus, List<Integer>>> activeCFNodes,
                                Map<Locus, Set<Conversion>> activeConversions) {

        for (Locus locus : acg.getConvertibleLoci()) {
            List<Integer> startSites = new ArrayList<>();
            List<Integer> endSites = new ArrayList<>();
            for (Node node : activeCFNodes.keySet()) {
                for (int i = 0; i < activeCFNodes.get(node).get(locus).size(); i += 2) {
                    startSites.add(activeCFNodes.get(node).get(locus).get(i));
                    endSites.add(activeCFNodes.get(node).get(locus).get(i + 1));
                }
            }

            for (Conversion conv : activeConversions.get(locus)) {
                for (int i = 0; i < affectedSites.get(conv).size(); i += 2) {
                    startSites.add(affectedSites.get(conv).get(i));
                    endSites.add(affectedSites.get(conv).get(i + 1));
                }
            }

            Collections.sort(startSites);
            Collections.sort(endSites);

            for (int i = 0; i < startSites.size() - 1; i++) {
                if (startSites.get(i + 1) < endSites.get(i))
                    return false;
            }
        }

        return true;
    }



}
