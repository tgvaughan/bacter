/*
 * Copyright (C) 2014 Tim Vaughan <tgvaughan@gmail.com>
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

package argbeast.util;

import argbeast.Recombination;
import argbeast.RecombinationGraph;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import com.google.common.collect.Lists;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/**
 * A home for methods of use mostly in debugging and unit testing.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class UtilMethods {
    
    /**
     * Create Alignment object representing alignment of a region
     * corresponding to a single marginal tree.
     * 
     * @param alignment
     * @param arg
     * @param recomb
     * @return
     * @throws Exception 
     */
    public static Alignment createMarginalAlignment(Alignment alignment,
            RecombinationGraph arg, Recombination recomb) throws Exception {
        List<Sequence> sequences = Lists.newArrayList();

        for (int leafIdx=0; leafIdx<alignment.getNrTaxa(); leafIdx++) {
            List<Integer> stateSequence;
            
            if (recomb == null) {
                stateSequence = Lists.newArrayList();
                
                // Portions of CF sequence before each converted region
                int i=0;
                for (int r=0; r<arg.getNRecombs(); r++) {
                    while (i<arg.getRecombinations().get(r+1).getStartLocus()) {
                        stateSequence.add(alignment.getCounts().get(leafIdx).get(i));
                        i += 1;
                    }
                    i=arg.getRecombinations().get(r+1).getEndLocus()+1;
                }
                
                // Any remaining CF sequence
                while (i<arg.getSequenceLength()) {
                    stateSequence.add(alignment.getCounts().get(leafIdx).get(i));
                    i += 1;
                }

            } else {

                stateSequence = alignment.getCounts().get(leafIdx)
                        .subList(recomb.getStartLocus(), recomb.getEndLocus()+1);
            }
            
            sequences.add(new Sequence(alignment.getTaxaNames().get(leafIdx),
                    alignment.getDataType().state2string(stateSequence)));
        }
        
        Alignment margAlignment = new Alignment(sequences,
                alignment.getDataType().getStateCount(),
                alignment.getDataType().getDescription());
        
        return margAlignment;
    }
    
    /**
     * Tests whether treeA and treeB are equivalent.  That is, whether they
     * have the same node heights (to within tolerance) and clade sets.
     * 
     * @param treeA
     * @param treeB
     * @param tolerance
     * @return 
     */
    public static boolean treesEquivalent(Tree treeA, Tree treeB, double tolerance) {

        // Early exit if trees are different sizes
        if (treeA.getLeafNodeCount() != treeB.getLeafNodeCount())
            return false;
        
        Map<Clade, Double> cladeHeightsA = getCladeHeights(treeA);
        Map<Clade, Double> cladeHeightsB = getCladeHeights(treeB);
        
        if (!cladeHeightsA.keySet().containsAll(cladeHeightsB.keySet()))
            return false;
        
        for (Clade clade : cladeHeightsA.keySet()) {
            if (Math.abs(cladeHeightsA.get(clade)-cladeHeightsB.get(clade))>tolerance)
                return false;
        }
        
        return true;
    }
    
    /**
     * Method to test whether the topologies of two trees are equivalent.
     * 
     * @param treeA
     * @param treeB
     * @return true iff topologies are equivalent.
     */
    public static boolean topologiesEquivalent(Tree treeA, Tree treeB) {
        
        // Early exit if trees are different sizes
        if (treeA.getLeafNodeCount() != treeB.getLeafNodeCount())
            return false;
        
        Map<Clade, Double> cladeHeightsA = getCladeHeights(treeA);
        Map<Clade, Double> cladeHeightsB = getCladeHeights(treeB);
        
        return cladeHeightsA.keySet().containsAll(cladeHeightsB.keySet());
    }
    
    /**
     * Retrieve clades and heights of clade MRCAs from tree.
     * 
     * @param tree
     * @return Map from clades to corresponding MRCA heights.
     */
    private static Map<Clade, Double> getCladeHeights(Tree tree) {
        Map<Clade, Double> cladeHeights = new HashMap<Clade, Double>();
        
        for (Node node : tree.getInternalNodes())
            cladeHeights.put(new Clade(node), node.getHeight());
        
        return cladeHeights;
    }
    
    /**
     * Convenience clade class.
     */
    private static class Clade extends HashSet<Integer> {

        /**
         * Construct clade from leaves below node.
         * 
         * @param node 
         */
        public Clade(Node node) {
            for (Node leaf : node.getAllLeafNodes())
                add(leaf.getNr());
        }
    };
}
