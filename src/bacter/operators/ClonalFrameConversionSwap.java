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

package bacter.operators;

import bacter.Conversion;
import beast.core.Description;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import java.util.Collections;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which swaps the role of the clonal frame with that of "
        + "one of the marginal trees resulting from a conversion.  Note that "
        + "this move conserves the marginal trees themselves.")
public class ClonalFrameConversionSwap extends EdgeCreationOperator {

    @Override
    public double proposal() {
        double logHR = 0.0;

        if (acg.getConvCount()==0 || getGapCount()==0)
            return Double.NEGATIVE_INFINITY;
        
        // Choose recombination
        int ridx = Randomizer.nextInt(acg.getConvCount())+1;
        Conversion recomb = acg.getConversions().get(ridx);

        if (recomb.getNode1()==recomb.getNode2())
            return Double.NEGATIVE_INFINITY;

        // Calculate probability of reverse move:
        logHR += getReverseMoveProb(recomb);
        
        // Perform swap and incorporate probability of forward move:
        logHR -= performSwap(recomb);
        
        return logHR;
    }
    
    /**
     * Return number of contiguous regions corresponding to clonal frame.
     * 
     * @return clonal frame region count
     */
    private int getGapCount() {
        int count = 0;

        if (acg.getConvCount()>0 && acg.getConversions().get(1).getStartSite()>0)
            count += 1;
        
        if (acg.getConvCount()>1)
            count += acg.getConvCount()-1;
        
        if (acg.getConvCount() == 0
                || acg.getConversions().get(acg.getConvCount()).getEndSite()<acg.getSequenceLength()-1)
            count += 1;
        
        return count;
    }
    
    /**
     * Perform swap operation.
     * 
     * @param recomb selected recombination
     * @return log probability of new state
     */
    protected double performSwap(Conversion recomb) {
        double logP = 0.0;
        
        // Choose unconverted (clonal frame) fragment:
        int gapCount = getGapCount();
        int chosenGapIdx = Randomizer.nextInt(gapCount);
        
        // Record details required to effect topology change:
        Node pivot = recomb.getNode1();
        Node floating = pivot.getParent();
        Node sister = floating.getLeft() == pivot
                ? floating.getRight()
                : floating.getLeft();
        
        // Create recombination corresponding to original clonal frame:
        Conversion oldCFrecomb = new Conversion();
        oldCFrecomb.setNode1(pivot);
        oldCFrecomb.setHeight1(recomb.getHeight1());
        if (recomb.getNode2()==sister)
            oldCFrecomb.setNode2(floating);
        else
            oldCFrecomb.setNode2(sister);
        oldCFrecomb.setHeight2(floating.getHeight());
        
        // Make marginal tree of chosen recomb the new clonal frame.

        floating.removeChild(sister);
        
        if (recomb.getNode2() == floating)
            recomb.setNode2(sister);
        
        if (floating.isRoot())
            sister.setParent(null);
        else {
            Node grandparent = floating.getParent();
            floating.setParent(null);
            grandparent.removeChild(floating);
            grandparent.addChild(sister);
        }
        
        Node newSister = recomb.getNode2();

        if (!newSister.isRoot()) {
            Node newParent = newSister.getParent();
            newParent.removeChild(newSister);
            newParent.addChild(floating);
        }
        floating.addChild(newSister);
        floating.setHeight(recomb.getHeight2());
        
        // Ensure any change of root is catered for
        
        if (floating.isRoot())
            acg.setRoot(floating);
        
        if (sister.isRoot())
            acg.setRoot(sister);
        
        // Record the site ranges the existing conversions apply to and
        // delete those conversions. (The CF now applies to those conversions.)
        
        List<Integer> startSites = Lists.newArrayList();
        List<Integer> endSites = Lists.newArrayList();

        while (acg.getConversions().size()>1) {
            Conversion thisRecomb = acg.getConversions().get(acg.getConvCount());
            startSites.add(thisRecomb.getStartSite());
            endSites.add(thisRecomb.getEndSite());
            acg.deleteConversion(acg.getConversions().get(acg.getConvCount()));
        }
        
        Collections.reverse(startSites);
        Collections.reverse(endSites);
        
        // Create new conversions corresponding to the regions originally
        // governed by the CF. One of these will be equivalent to the old CF.
        
        int gapIdx = 0;
        if (startSites.get(0)>0) {
            Conversion newRecomb;

            if (chosenGapIdx==0) {
                newRecomb = oldCFrecomb;
            } else {
                newRecomb = new Conversion();
                logP += attachEdge(newRecomb);                
            }
            
            newRecomb.setStartSite(0);
            newRecomb.setEndSite(startSites.get(0)-1);
            acg.addConversion(newRecomb);
            
            gapIdx += 1;
        }
        
        for (int i=0; i<startSites.size()-1; i++) {
            Conversion newRecomb;
            if (chosenGapIdx==gapIdx) {
                newRecomb = oldCFrecomb;
            } else {
                newRecomb = new Conversion();
                logP += attachEdge(newRecomb);
            }
            newRecomb.setStartSite(endSites.get(i)+1);
            newRecomb.setEndSite(startSites.get(i+1)-1);
            acg.addConversion(newRecomb);
            
            gapIdx += 1;
        }
        
        if (endSites.get(endSites.size()-1)<acg.getSequenceLength()-1) {
            Conversion newRecomb;
            if (chosenGapIdx==gapIdx) {
                newRecomb = oldCFrecomb;
            } else {
                newRecomb = new Conversion();
                logP += attachEdge(newRecomb);
            }
            newRecomb.setStartSite(endSites.get(endSites.size()-1)+1);
            newRecomb.setEndSite(acg.getSequenceLength()-1);
            acg.addConversion(newRecomb);
            
            gapIdx += 1;
        }
        
        return logP;
    }
    
    /**
     * Obtain probability of reverse move.
     * 
     * @param recomb recombination chosen for forward move
     * @return log of reverse move probability.
     */
    protected double getReverseMoveProb(Conversion recomb) {
        double logP = 0.0;

        // Probability of drawing existing recombinant edges in x from
        // the clonal frame in x'.  (This excludes the conversion corresponding
        // to the chosen gap in x'.)
        for (Conversion thisRecomb : acg.getConversions()) {
            if (thisRecomb == null || thisRecomb == recomb)
                continue;
            
            logP += getEdgeAttachmentProb(thisRecomb);
        }
        
        return logP;
    }
}
