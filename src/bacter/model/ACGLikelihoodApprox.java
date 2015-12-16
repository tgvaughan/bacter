package bacter.model;

import bacter.*;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Node;
import beast.math.Binomial;

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

    int nLeaves;
    int[][] cumulativeHD;
    int[] pairNrLookup;
    Alignment alignment;
    ConversionGraph acg;
    Locus locus;

    public ACGLikelihoodApprox() { }

    @Override
    public void initAndValidate() throws Exception {

        alignment = alignmentInput.get();
        acg = acgInput.get();
        locus = locusInput.get();

        computePairwiseDistances();
    }


    @Override
    public double calculateLogP() throws Exception {
        logP = 0.0;

        Map<Double, Coalescence> heightMap = getCoalescenceHeights();

        for (Map.Entry<Double, Coalescence> entry : heightMap.entrySet()) {
            double height = entry.getKey();
            Coalescence coalescence = entry.getValue();

            for (int i=0; i<coalescence.getIntervalCount(); i++) {
                BitSet dl1 = coalescence.descendantLeaves1.get(i);
                BitSet dl2 = coalescence.descendantLeaves2.get(i);
                int x = coalescence.siteRanges.get(2*i);
                int y = coalescence.siteRanges.get(2*i + 1);

                for (int nr1 = dl1.nextSetBit(0); nr1>=0; nr1 = dl1.nextSetBit(nr1+1)) {
                    for (int nr2 = dl2.nextSetBit(0); nr2>=0; nr2 = dl2.nextSetBit(nr2+1)) {
                        double time = 2*height
                                - acg.getNode(nr1).getHeight()
                                - acg.getNode(nr2).getHeight();

                        int h = getPairwiseDistance(nr1, nr2, x, y);

                        logP += getHDProbability(h, time, y-x);
                    }
                }

            }
        }

        return logP;
    }

    /**
     * Returns the (log) probability of observing h segregating sites out of
     * a total of siteCount sites when the sequence is left to evolve for
     * a given time.
     *
     * @param h observed Hamming distance
     * @param time total time for evolution
     * @param siteCount total number of sites (h<=siteCount)
     * @return log probability
     */
    public double getHDProbability(int h, double time, int siteCount) {
        double p = 0.75*(1-Math.exp(-4.0/3.0*time*substRateInput.get().getValue()));
        return Binomial.logChoose(siteCount, h)
                + h*Math.log(p) + (siteCount-h)*Math.log(1.0-p);
    }

    /**
     * @return map from heights of coalescences to objects describing
     * the sites and samples they involve.
     */
    public Map<Double, Coalescence> getCoalescenceHeights() {

        Map<Double, Coalescence> heightMap = new HashMap<>();

        Map<Node, SiteAncestry> activeCFNodes = new HashMap<>();
        Map<Conversion, SiteAncestry> activeConversions = new HashMap<>();

        ACGEventList acgEventList = new ACGEventList(acg, locus);

        for (ACGEventList.Event event : acgEventList.getACGEvents()) {

            switch (event.type) {
                case CF_LEAF:
                    activeCFNodes.put(event.node, new SiteAncestry(event.node, locus));

                    break;

                case CF_COALESCENCE:
                    Node node1 = event.node.getLeft();
                    Node node2 = event.node.getRight();

                    SiteAncestry ancestryCF = new SiteAncestry();
                    Coalescence coalescenceCF = new Coalescence();
                    activeCFNodes.get(node1).merge(activeCFNodes.get(node2),
                            coalescenceCF, ancestryCF);

                    activeCFNodes.remove(node1);
                    activeCFNodes.remove(node2);
                    activeCFNodes.put(event.node, ancestryCF);

                    if (coalescenceCF.getIntervalCount()>0)
                        heightMap.put(event.t, coalescenceCF);

                    break;

                case CONV_DEPART:
                    SiteAncestry inside = new SiteAncestry();
                    SiteAncestry outside = new SiteAncestry();
                    activeCFNodes.get(event.node).split(
                            event.conversion.getStartSite(),
                            event.conversion.getEndSite()+1,
                            inside, outside);

                    activeCFNodes.put(event.node, outside);
                    activeConversions.put(event.conversion, inside);

                    break;

                case CONV_ARRIVE:

                    SiteAncestry ancestry = new SiteAncestry();
                    Coalescence coalescence = new Coalescence();
                    activeCFNodes.get(event.node).merge(activeConversions.get(event.conversion),
                            coalescence, ancestry);

                    activeCFNodes.put(event.node, ancestry);
                    activeConversions.remove(event.conversion);

                    if (coalescence.getIntervalCount()>0)
                        heightMap.put(event.t, coalescence);

                    break;
            }

        }

        return heightMap;
    }

    public void computePairwiseDistances() {
        // Pre-compute pairwise distance tables

        nLeaves = acg.getLeafNodeCount();
        int nPairs = nLeaves*(nLeaves-1)/2;
        cumulativeHD = new int[nPairs][alignment.getSiteCount()+1];
        pairNrLookup = new int[nLeaves*nLeaves];

        for (int siteBoundary=0; siteBoundary<alignment.getSiteCount()+1; siteBoundary++) {
            int pair = 0;
            for (int tIdx1=0; tIdx1<nLeaves; tIdx1++) {
                for (int tIdx2=tIdx1+1; tIdx2<nLeaves; tIdx2++) {
                    pairNrLookup[tIdx1 + tIdx2*nLeaves] = pair;
                    pairNrLookup[tIdx2 + tIdx1*nLeaves] = pair;
                    if (siteBoundary == 0)
                        cumulativeHD[pair][siteBoundary] = 0;
                    else
                        cumulativeHD[pair][siteBoundary] = cumulativeHD[pair][siteBoundary - 1];

                    if (siteBoundary > 0) {
                        int patternIdx = alignment.getPatternIndex(siteBoundary - 1);

                        if (alignment.getPattern(tIdx1, patternIdx)
                                != alignment.getPattern(tIdx2, patternIdx))
                            cumulativeHD[pair][siteBoundary] += 1;
                    }

                    pair += 1;
                }
            }
        }
    }

    public int getPairwiseDistance(int node1Nr, int node2Nr, int x, int y) {
        int pairNr = pairNrLookup[node1Nr*nLeaves + node2Nr];
        return cumulativeHD[pairNr][y] - cumulativeHD[pairNr][x];
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
