package bacter.model;

import bacter.*;
import bacter.util.IntRanges;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Node;

import java.util.*;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Approximate ACG likelihood using pairwise distances.")
public class ACGLikelihoodApprox extends Distribution {

    public Input<ConversionGraph> acgInput = new Input<>(
            "acg",
            "Conversion graph.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> substRateInput = new Input<>(
            "substitutionRate",
            "Substitution rate.",
            Input.Validate.REQUIRED);

    public Input<Alignment> alignmentInput = new Input<>(
            "alignment",
            "Alignment corresponding to particular locus.",
            Input.Validate.REQUIRED);

    public Input<Locus> locusInput = new Input<>(
            "locus",
            "Locus alignment is associated with.",
            Input.Validate.REQUIRED);

    int nPairs;
    int[][] cumulativeHD;
    Alignment alignment;
    ConversionGraph acg;
    Locus locus;

    public ACGLikelihoodApprox() { }

    @Override
    public void initAndValidate() throws Exception {

        alignment = alignmentInput.get();
        acg = acgInput.get();
        locus = locusInput.get();

        // Pre-compute pairwise distance tables

        int nLeaves = acg.getLeafNodeCount();
        nPairs = nLeaves*(nLeaves+1)/2;
        cumulativeHD = new int[nPairs][alignment.getSiteCount()];

        for (int site=0; site<alignment.getSiteCount(); site++) {
            int pair = 0;
            for (int tIdx1=1; tIdx1<nLeaves; tIdx1++) {
                for (int tIdx2=0; tIdx2<tIdx1; tIdx2++) {
                    if (site==0)
                        cumulativeHD[pair][site] = 0;
                    else
                        cumulativeHD[pair][site] = cumulativeHD[pair][site-1];

                    int patternIdx = alignment.getPatternIndex(site);

                    if (alignment.getPattern(tIdx1, patternIdx)
                            != alignment.getPattern(tIdx2, patternIdx))
                        cumulativeHD[pair][site] += 1;

                    pair += 1;
                }
            }
        }
    }



    @Override
    public double calculateLogP() throws Exception {
        logP = 0.0;

        return logP;
    }

    Map<Double, List<Integer>> getCoalescenceHeights() {
        Map<Double, List<Integer>> heightMap = new HashMap<>();

        Map<Node, List<Integer>> activeCFNodes = new HashMap<>();
        Map<Conversion, List<Integer>> activeConversions = new HashMap<>();

        ACGEventList acgEventList = new ACGEventList(acg, locus);

        for (ACGEventList.Event event : acgEventList.getACGEvents()) {

            switch (event.type) {
                case CF_LEAF:
                    List<Integer> leafSites = new ArrayList<>();
                    leafSites.add(0);
                    leafSites.add(locus.getSiteCount()-1);
                    activeCFNodes.put(event.node, leafSites);

                    break;

                case CF_COALESCENCE:
                    Node node1 = event.node.getLeft();
                    Node node2 = event.node.getRight();

                    List<Integer> ancestralSitesCF =
                            IntRanges.getUnion(activeCFNodes.get(node1),
                                        activeCFNodes.get(node2));

                    List<Integer> coalescingSitesCF =
                            IntRanges.getIntersection(
                                    activeCFNodes.get(node1),
                                    activeCFNodes.get(node2));

                    if (!coalescingSitesCF.isEmpty())
                        heightMap.put(event.node.getHeight(), coalescingSitesCF);

                    activeCFNodes.remove(node1);
                    activeCFNodes.remove(node2);
                    activeCFNodes.put(event.node, ancestralSitesCF);

                    break;

                case CONV_DEPART:
                    List<Integer> inside = new ArrayList<>();
                    List<Integer> outside = new ArrayList<>();
                    IntRanges.partitionRanges(activeCFNodes.get(event.node),
                            event.conversion.getStartSite(),
                            event.conversion.getEndSite()+1,
                            inside, outside);

                    activeCFNodes.put(event.node, outside);
                    activeConversions.put(event.conversion, outside);

                    break;

                case CONV_ARRIVE:
                    activeCFNodes.put(event.node,
                            IntRanges.getUnion(
                                    activeConversions.get(event.conversion),
                                    activeCFNodes.get(event.node)));

                    List<Integer> coalescingSites =
                            IntRanges.getIntersection(
                                    activeConversions.get(event.conversion),
                                    activeCFNodes.get(event.node));

                    if (!coalescingSites.isEmpty())
                        heightMap.put(event.conversion.getHeight2(), coalescingSites);

                    activeConversions.remove(event.conversion);

                    break;
            }

        }

        return heightMap;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Sampling from this " +
                "ACGLikelihoodApprox is not supported.");
    }
}
