/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bacter.acgannotator;

import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.Locus;
import beast.base.evolution.tree.Node;
import beastfx.app.treeannotator.CladeSystem;
import beast.base.util.Randomizer;

import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

/**
 * Adds conversion summary tools to CladeSystem.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ACGCladeSystem extends CladeSystem {

    protected Map<BitSetPair, Map<Locus, List<Conversion>>> conversionLists = new HashMap<>();
    protected Map<BitSetPair, List<Conversion>> conversionListsTemp = new HashMap<>();
    protected List<Map<BitSet, Map<BitSet, Long>>> geneFlow = new ArrayList<>();
    protected BitSet[] bitSets;

    protected int acgIndex = 1;

    public ACGCladeSystem() { }

    public ACGCladeSystem(ConversionGraph acg) {
        add(acg, true);
    }

    /**
     * Assemble list of bitSets for this ACG.
     */
    public BitSet[] getBitSets(ConversionGraph acg) {

        if (bitSets == null)
            bitSets = new BitSet[acg.getNodeCount()];

        applyToClades(acg.getRoot(), (cladeNode, bits) -> {
            bitSets[cladeNode.getNr()] = bits;
            return null;
        });

        return bitSets;
    }

    /**
     * Add conversions described on provided acg to the internal list
     * for later summary.
     *
     * @param acg conversion graph from which to extract conversions
     * @param receiverBranchMode
     */
    public void collectConversions(ConversionGraph acg, boolean receiverBranchMode) {
        getBitSets(acg);

        Map<BitSet,Map<BitSet,Long>> geneFlowTemp = new HashMap<>();

        // Assemble list of conversions for each pair of clades on each locus
        //receiver branch mode edit
        for (Locus locus : acg.getConvertibleLoci()) {

            conversionListsTemp.clear();
            for (Conversion conv : acg.getConversions(locus))  {
                conv.acgIndex = acgIndex;
                BitSetPair bsPair = new BitSetPair(conv);
                if(receiverBranchMode)
                    bsPair.to = new BitSet(); //set the to bitSet to empty if the receiverBranch mode is used

                if (!conversionListsTemp.containsKey(bsPair))
                    conversionListsTemp.put(bsPair, new ArrayList<>());

                conversionListsTemp.get(bsPair).add(conv);

                // Record gene flow
                if (!geneFlowTemp.containsKey(bsPair.from))
                    geneFlowTemp.put(bsPair.from, new HashMap<>());

                long oldFlow = 0;
                if (geneFlowTemp.get(bsPair.from).containsKey(bsPair.to))
                    oldFlow = geneFlowTemp.get(bsPair.from).get(bsPair.to);

                geneFlowTemp.get(bsPair.from).put(bsPair.to, oldFlow + conv.getSiteCount());
            }

            // Merge overlapping conversions:
            for (BitSetPair bsPair : conversionListsTemp.keySet()) {
                List<Conversion> merged = mergeOverlappingConvs(
                            conversionListsTemp.get(bsPair));

                if (!conversionLists.containsKey(bsPair))
                    conversionLists.put(bsPair, new HashMap<>());
                if (!conversionLists.get(bsPair).containsKey(locus))
                    conversionLists.get(bsPair).put(locus, new ArrayList<>());

                conversionLists.get(bsPair).get(locus).addAll(merged);

            }
        }

        geneFlow.add(geneFlowTemp);

        acgIndex += 1;
    }

    private List<Conversion> mergeOverlappingConvs(List<Conversion> conversions) {
        List<Conversion> mergedList = new ArrayList<>();

        List<Conversion> convOrderedByStart = new ArrayList<>(conversions);
        convOrderedByStart.sort((o1, o2) -> o1.getStartSite() - o2.getStartSite());

        List<Conversion> convOrderedByEnd = new ArrayList<>(conversions);
        convOrderedByEnd.sort((o1, o2) -> o1.getEndSite() - o2.getEndSite());


        int nActive = 0;
        Conversion currentMergedConv = null;
        int mergedConvCount = 0;
        List<Double> mergedConvHeight1 = new ArrayList<>();
        List<Double> mergedConvHeight2 = new ArrayList<>();
        List<Double> mergedConvHeight1First = new ArrayList<Double>();
        List<Double> mergedConvHeight2First = new ArrayList<Double>();
        List<Node> mergedNodes2 = new ArrayList<>();
        List<Node> mergedNodes2First = new ArrayList<>();

        //circular genome mode edit
        Map<Node, Integer> node2counts = new HashMap<>();
        Map<Node, Integer> node2countsFirst = new HashMap<>();
        int numFirst = 0;

        boolean firstStep = true;
        int minOverlapStart = Integer.MAX_VALUE;
        int indConv = 0;
        for (int i = 0; i < convOrderedByEnd.size(); i++) {
            Conversion conv = convOrderedByStart.get(indConv);
            if (conv.getEndSite() < conv.getStartSite()) {
                nActive += 1;
                mergedConvCount += 1;
                mergedConvHeight1.add(conv.getHeight1());
                mergedConvHeight2.add(conv.getHeight2());
                mergedNodes2.add(conv.getNode2());
                minOverlapStart = Math.min(minOverlapStart, conv.getStartSite());
                currentMergedConv = conv.getStartSite() <= minOverlapStart ? conv.getCopy() : currentMergedConv;
                currentMergedConv.acgIndex = conv.acgIndex;
                convOrderedByStart.remove(indConv);
                indConv -= 1;
            }
            indConv += 1;
        }

        while (!convOrderedByStart.isEmpty() || !convOrderedByEnd.isEmpty()) {

            int nextStart = convOrderedByStart.isEmpty()
                    ? Integer.MAX_VALUE
                    : convOrderedByStart.get(0).getStartSite();

            int nextEnd = convOrderedByEnd.isEmpty()
                    ? Integer.MAX_VALUE
                    : convOrderedByEnd.get(0).getEndSite();

            if (nextStart < nextEnd) {
                nActive += 1;

                if (nActive == 1) {
                    currentMergedConv = convOrderedByStart.get(0).getCopy();
                    currentMergedConv.acgIndex = convOrderedByStart.get(0).acgIndex;
                    mergedConvCount = 1;
                    mergedConvHeight1.clear();
                    mergedConvHeight2.clear();
                    mergedNodes2.clear();
                    mergedConvHeight1.add(currentMergedConv.getHeight1());
                    mergedConvHeight2.add(currentMergedConv.getHeight2());
                    mergedNodes2.add(currentMergedConv.getNode2());
                } else {
                    mergedConvCount += 1;
                    mergedConvHeight1.add(convOrderedByStart.get(0).getHeight1());
                    mergedConvHeight2.add(convOrderedByStart.get(0).getHeight2());
                    mergedNodes2.add(convOrderedByStart.get(0).getNode2());
                }

                convOrderedByStart.remove(0);

            } else {
                nActive -= 1;

                if (nActive == 0 ) {
                    assert currentMergedConv != null;
                    currentMergedConv.setEndSite(nextEnd);
                    //receiver branch mode edit
                    node2counts.clear();
                    for (Node node2 : mergedNodes2) {
                        Integer count = node2counts.get(node2);
                        node2counts.put(node2, count != null ? count+1 : 1);
                    }
                    int maxNode2count = 0;
                    Node selectedNode = null;
                    for (Node node2 : node2counts.keySet()){
                        if (node2counts.get(node2) > maxNode2count){
                            maxNode2count = node2counts.get(node2);
                            selectedNode = node2;
                        } else if (node2counts.get(node2) == maxNode2count){
                            selectedNode = Randomizer.nextBoolean() ? node2 : selectedNode;
                        }
                    }
                    double sumSelectedHeights1 = 0;
                    double sumSelectedHeights2 = 0;
                    for (int i = 0; i <  mergedConvCount; i++ ){
                        if (mergedNodes2.get(i).equals(selectedNode))
                            sumSelectedHeights2 += mergedConvHeight2.get(i);
                        sumSelectedHeights1 += mergedConvHeight1.get(i);
                    }
                    currentMergedConv.setHeight1(sumSelectedHeights1 / mergedConvCount);
                    currentMergedConv.setHeight2(sumSelectedHeights2 / maxNode2count);
                    currentMergedConv.setNode2(selectedNode);
                    mergedList.add(currentMergedConv);
                    //circular genome mode edit
                    if (firstStep) {
                        mergedConvHeight1First = new ArrayList<Double>(mergedConvHeight1);
                        mergedConvHeight2First = new ArrayList<Double>(mergedConvHeight2);
                        mergedNodes2First = new ArrayList<Node>(mergedNodes2);
                        node2counts.forEach((key, value) -> node2countsFirst.merge(key, value, Integer::sum));
                        numFirst = mergedConvCount;
                        firstStep = false;
                    }
                }

                convOrderedByEnd.remove(0);
            }
        }
        //circular genome mode edit
        if (mergedList.size() > 1 && (currentMergedConv.getEndSite() >= minOverlapStart)) {
            mergedList.remove(mergedList.size()-1);
            mergedList.get(0).setStartSite(currentMergedConv.getStartSite());
            node2countsFirst.forEach((key, value) -> node2counts.merge(key, value, Integer::sum));
            mergedConvHeight1.addAll(mergedConvHeight1First);
            mergedConvHeight2.addAll(mergedConvHeight2First);
            mergedNodes2.addAll(mergedNodes2First);
            //receiver branch mode edit
            int maxNode2count = 0;
            Node selectedNode = null;
            for (Node node2 : node2counts.keySet()){
                if (node2counts.get(node2) > maxNode2count){
                    maxNode2count = node2counts.get(node2);
                    selectedNode = node2;
                } else if (node2counts.get(node2) == maxNode2count){
                    selectedNode = Randomizer.nextBoolean() ? node2 : selectedNode;
                }
            }
            double sumSelectedHeights1 = 0;
            double sumSelectedHeights2 = 0;
            for (int i = 0; i <  mergedConvCount; i++ ){
                if (mergedNodes2.get(i).equals(selectedNode))
                    sumSelectedHeights2 += mergedConvHeight2.get(i);
                sumSelectedHeights1 += mergedConvHeight1.get(i);
            }
            mergedList.get(0).setHeight1(sumSelectedHeights1 / (mergedConvCount + numFirst));
            mergedList.get(0).setHeight2(sumSelectedHeights2 / (maxNode2count));
            mergedList.get(0).setNode2(selectedNode);
        }
        return mergedList;
    }

    /**
     * Determine contiguous regions on specified locus where the fraction of
     * ACGs having a conversion active is greater than the given threshold.
     *
     * @param from BitSet representing source clade
     * @param to BitSet representing destination clade
     * @param locus locus to consider
     * @param threshold minimum fraction of sampled conversions included
     * @return List of regions
     */
    public List<ConversionSummary> getConversionSummaries(BitSet from, BitSet to,
                                                          Locus locus,
                                                          int nACGs,
                                                          double threshold) {

        BitSetPair bsPair = new BitSetPair(from, to);

        List<ConversionSummary> convSummaryList = new ArrayList<>();

        // Return empty list if no conversions meet the criteria.
        if (!conversionLists.containsKey(bsPair)
                || !conversionLists.get(bsPair).containsKey(locus))
            return convSummaryList;

        int thresholdCount = (int)Math.ceil(nACGs*threshold);

        List<Conversion> convOrderedByStart = new ArrayList<>();
        convOrderedByStart.addAll(conversionLists.get(bsPair).get(locus));
        convOrderedByStart.sort((Conversion o1, Conversion o2) ->
                o1.getStartSite() - o2.getStartSite());

        List<Conversion> convOrderedByEnd = new ArrayList<>();
        convOrderedByEnd.addAll(conversionLists.get(bsPair).get(locus));
        convOrderedByEnd.sort((Conversion o1, Conversion o2) ->
                o1.getEndSite() - o2.getEndSite());

        List<Conversion> activeConversions = new ArrayList<>();
        ConversionSummary conversionSummary = null;

        BitSet includedACGindices = new BitSet();

        //circular genome mode edit
        int numOverlap = 0;
        for (Conversion conv : convOrderedByStart) {
            if (conv.getEndSite() < conv.getStartSite()) {
                activeConversions.add(conv);
                includedACGindices.set(conv.acgIndex);
                numOverlap += 1;
            }
        }

        int overlapStartBound = Integer.MAX_VALUE;
        boolean overlapRegion = false;
        if (!activeConversions.isEmpty() && activeConversions.size() >= thresholdCount) {
            overlapStartBound = activeConversions.get(thresholdCount > 0 ? thresholdCount - 1 : 0).getStartSite();
            overlapRegion = true;
            conversionSummary = new ConversionSummary();
            convSummaryList.add(conversionSummary);
            conversionSummary.addConvs(activeConversions);
        }
        //circular genome mode edit
        int maxEndSite = !convOrderedByEnd.isEmpty() ? convOrderedByEnd.get(convOrderedByEnd.size() - 1).getEndSite() : Integer.MAX_VALUE;

        while (!convOrderedByStart.isEmpty() || !convOrderedByEnd.isEmpty()) {

            int nextStart = convOrderedByStart.isEmpty()
                    ? Integer.MAX_VALUE
                    : convOrderedByStart.get(0).getStartSite();

            int nextEnd = convOrderedByEnd.isEmpty()
                    ? Integer.MAX_VALUE
                    : convOrderedByEnd.get(0).getEndSite();

            if (nextStart < nextEnd) {
                activeConversions.add(convOrderedByStart.get(0));

                if (activeConversions.size() >= thresholdCount) {
                    if ( conversionSummary == null) {
                        conversionSummary = new ConversionSummary();
                        convSummaryList.add(conversionSummary);
                        conversionSummary.addConvs(activeConversions);

                        includedACGindices.clear();
                        for (Conversion conv : activeConversions)
                            includedACGindices.set(conv.acgIndex);
                    } else {
                        if (convOrderedByStart.get(0).getStartSite() <= convOrderedByStart.get(0).getEndSite()) {
                            conversionSummary.addConv(convOrderedByStart.get(0));
                            includedACGindices.set(convOrderedByStart.get(0).acgIndex);
                        }
                    }
                }
                convOrderedByStart.remove(0);
            } else {
                //circular genome mode edit
                if (conversionSummary != null && overlapRegion && nextEnd > overlapStartBound && nextEnd == maxEndSite) {
                    if (convSummaryList.size() > 1) {
                        convSummaryList.remove(conversionSummary);
                        convSummaryList.get(0).mergeConvSum(conversionSummary);
                        convSummaryList.get(0).nIncludedACGs += conversionSummary.nIncludedACGs;
                    }
                    conversionSummary = null;
                }
                activeConversions.remove(convOrderedByEnd.get(0));
                if (activeConversions.size() == thresholdCount-1) {
                    assert conversionSummary != null;
                    conversionSummary.nIncludedACGs = includedACGindices.cardinality();
                    conversionSummary = null;
                }
                convOrderedByEnd.remove(0);
            }
        }
        if (conversionSummary != null && thresholdCount > 0) {
            convSummaryList.remove(conversionSummary);
        }
        return convSummaryList;
    }

    /**
     * @return list of maps specifying gene flow between clades.
     */
    public List<Map<BitSet,Map<BitSet,Long>>> getGeneFlowMaps() {
        return geneFlow;
    }

    /**
     * Apply a function to each sub-clade.
     *
     * @param node MRCA of clade
     * @param function function to apply. Given sub-clade parent node
     *                 and bitset as arguments.
     * @return BitSet representing clade.
     */
    public BitSet applyToClades(Node node, BiFunction<Node, BitSet, Void> function) {
        BitSet bits = new BitSet();

        if (node.isLeaf()) {
            bits.set(2 * getTaxonIndex(node));
        } else {
            for (Node child : node.getChildren())
                bits.or(applyToClades(child, function));
        }

        function.apply(node, bits);

        return bits;
    }


        /**
         * Get bitset corresponding to a clade
         *
         * @param node MRCA of clade
         * @return BitSet representing clade.
         */
        //receiver branch mode edit
        public BitSet getBitSet(Node node) {
            BitSet bitset = new BitSet();

            if (node.isLeaf()) {
                bitset.set(2 * node.getNr());
            } else {
                for (Node child : node.getChildren())
                    bitset.or(getBitSet(child));
            }

            return bitset;
        }



    /**
     * Class representing an ordered pair of BitSets.
     */
    protected class BitSetPair {
        public BitSet from, to;

        public BitSetPair(BitSet from, BitSet to) {
            this.from = from;
            this.to = to;
        }

        public BitSetPair(Conversion conv) {
            this.from = bitSets[conv.getNode1().getNr()];
            this.to = bitSets[conv.getNode2().getNr()];
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            BitSetPair that = (BitSetPair) o;

            return from.equals(that.from) && to.equals(that.to);

        }

        @Override
        public int hashCode() {
            int result = from.hashCode();
            result = 31 * result + to.hashCode();
            return result;
        }

        @Override
        public String toString() {
            return from.toString() + " -> " + to.toString();
        }
    }

    /**
     * Class representing a summary of similar conversions between two
     * points in the summarized clonal frame.
     */
    public class ConversionSummary {

        List<Double> height1s = new ArrayList<>();
        List<Double> height2s = new ArrayList<>();
        List<Integer> startSites = new ArrayList<>();
        List<Integer> ends = new ArrayList<>();
        //receiver branch mode edit
        List<BitSet> node2s = new ArrayList<>();
        public int nIncludedACGs = 0;

        /**
         * Add metrics associated with given conversion to summary.
         *
         * @param conv conversion
         */
        public void addConv(Conversion conv) {
            height1s.add(conv.getHeight1());
            height2s.add(conv.getHeight2());
            startSites.add(conv.getStartSite());
            ends.add(conv.getEndSite());
            node2s.add(getBitSet(conv.getNode2()));
        }

        /**
         * Add metrics associated with each of the conversions in the
         * given list to the summary.
         *
         * @param convs list of conversions
         */
        public void addConvs(List<Conversion> convs) {
            for (Conversion conv : convs)
                addConv(conv);
        }

        /**
         * Merge metrics associated in given conversion summary with this summary.
         *
         //* @param ConversionSummary convSum
         */
        public void mergeConvSum(ConversionSummary convSum) {
            height1s.addAll(convSum.height1s);
            height2s.addAll(convSum.height2s);
            startSites.addAll(convSum.startSites);
            ends.addAll(convSum.ends);
            node2s.addAll(convSum.node2s);
        }

        /**
         * @return number of conversions included in summary.
         */
        public int summarizedConvCount() {
            return height1s.size();
        }
    }
}
