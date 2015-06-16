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

package bacter.util;

import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.Locus;
import bacter.Region;
import beast.app.treeannotator.CladeSystem;
import beast.evolution.tree.Node;

import java.util.*;
import java.util.function.BiFunction;

/**
 * Adds conversion summary tools to CladeSystem.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ACGCladeSystem extends CladeSystem {

    protected Map<BitSetPair, Map<Locus, List<Conversion>>> conversionLists = new HashMap<>();
    protected BitSet[] bitSets;

    public ACGCladeSystem() { }

    public ACGCladeSystem(ConversionGraph acg) {
        add(acg, true);
    }

    /**
     * Assemble list of bitSets for this ACG.
     */
    protected void getBitSets(ConversionGraph acg) {

        if (bitSets == null)
            bitSets = new BitSet[acg.getNodeCount()];

        applyToClades(acg.getRoot(), (cladeNode, bits) -> {
            bitSets[cladeNode.getNr()] = bits;
            return null;
        });

    }

    /**
     * Add conversions described on provided acg to the internal list
     * for later summary.
     *
     * @param acg conversion graph from which to extract conversions
     */
    public void collectConversions(ConversionGraph acg) {
        getBitSets(acg);

        conversionLists.clear();

        // Assemble list of conversions for each pair of clades on each locus
        for (Locus locus : acg.getLoci()) {
            for (Conversion conv : acg.getConversions(locus))  {
                BitSetPair bsPair = new BitSetPair(conv);

                if (!conversionLists.containsKey(bsPair))
                    conversionLists.put(bsPair, new HashMap<>());
                if (!conversionLists.get(bsPair).containsKey(locus))
                    conversionLists.get(bsPair).put(locus, new ArrayList<>());

                conversionLists.get(bsPair).get(locus).add(conv);
            }

            // Merge overlapping conversions:
            for (BitSetPair bsPair : conversionLists.keySet()) {
                if (conversionLists.get(bsPair).containsKey(locus)) {
                    List<Conversion> merged = mergeOverlappingConvs(
                            conversionLists.get(bsPair).get(locus));

                    conversionLists.get(bsPair).put(locus, merged);
                }
            }
        }
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
        double mergedConvHeight1 = 0.0;
        double mergedConvHeight2 = 0.0;

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
                    mergedConvCount = 1;
                    mergedConvHeight1 = currentMergedConv.getHeight1();
                    mergedConvHeight2 = currentMergedConv.getHeight2();
                } else {
                    mergedConvCount += 1;
                    mergedConvHeight1 += convOrderedByStart.get(0).getHeight1();
                    mergedConvHeight2 += convOrderedByStart.get(0).getHeight2();
                }

                convOrderedByStart.remove(0);

            } else {
                nActive -= 1;

                if (nActive == 0 ) {
                    assert currentMergedConv != null;
                    currentMergedConv.setEndSite(nextEnd);
                    currentMergedConv.setHeight1(mergedConvHeight1/mergedConvCount);
                    currentMergedConv.setHeight2(mergedConvHeight2 / mergedConvCount);
                    mergedList.add(currentMergedConv);
                }

                convOrderedByEnd.remove(0);
            }
        }

        return mergedList;
    }

    /**
     * Determine contiguous regions on specified locus where the fraction of
     * ACGs having a conversion active is greater than the given threshold.
     *
     * @param bsPair location of conversions on tree
     * @param locus locus to consider
     * @param threshold minimum fraction of sampled conversions included
     * @return List of regions
     */
    public List<Region> getSignificantRegions(BitSetPair bsPair, Locus locus, double threshold) {

        List<Region> regions = new ArrayList<>();

        List<Conversion> convOrderedByStart = new ArrayList<>();
        convOrderedByStart.addAll(conversionLists.get(bsPair).get(locus));
        convOrderedByStart.sort((Conversion o1, Conversion o2) ->
                o1.getStartSite() - o2.getStartSite());

        List<Conversion> convOrderedByEnd = new ArrayList<>();
        convOrderedByEnd.addAll(conversionLists.get(bsPair).get(locus));
        convOrderedByEnd.sort((Conversion o1, Conversion o2) ->
                o1.getEndSite() - o2.getEndSite());

        Set<Conversion> activeConversions = new HashSet<>();

        int lastBoundary = 0;

        while (!convOrderedByStart.isEmpty() || !convOrderedByEnd.isEmpty()) {

            int nextStart;
            if (!convOrderedByStart.isEmpty())
                nextStart = convOrderedByStart.get(0).getStartSite();
            else
                nextStart = Integer.MAX_VALUE;

            int nextEnd;
            if (!convOrderedByEnd.isEmpty())
                nextEnd = convOrderedByEnd.get(0).getEndSite() + 1;
            else
                nextEnd = Integer.MAX_VALUE;

            int nextBoundary = Math.min(nextStart, nextEnd);
            if (nextBoundary > lastBoundary) {
                Region region = new Region();
                region.leftBoundary = lastBoundary;
                region.rightBoundary = nextBoundary;
                region.activeConversions.addAll(activeConversions);
                regions.add(region);
            }

            if (nextStart < nextEnd) {
                activeConversions.add(convOrderedByStart.get(0));
                convOrderedByStart.remove(0);
                lastBoundary = nextStart;
            } else {
                activeConversions.remove(convOrderedByEnd.get(0));
                convOrderedByEnd.remove(0);
                lastBoundary = nextEnd;
            }
        }

        if (lastBoundary < locus.getSiteCount()) {
            Region region = new Region();
            region.leftBoundary = lastBoundary;
            region.rightBoundary = locus.getSiteCount();
            regions.add(region);
        }

        return null;
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
     * Class representing an ordered pair of BitSets.
     */
    public class BitSetPair {
        public BitSet from, to;

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
}
