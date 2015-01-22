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

package bacter.operators.restricted;

import bacter.operators.ConversionGraphOperator;
import bacter.Conversion;
import beast.core.Description;
import beast.core.Input;
import beast.util.Randomizer;
import feast.input.In;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which in one direction splits an existing recombination "
        + "in two and in the other merges two adjacent recombinations.")
public class MergeSplitConversion extends ConversionGraphOperator {
    
    public Input<Double> gapSizeInput = new In<Double>("expectedGapSize",
            "Expected gap size.").setDefault(10.0);
    
    private double gapRate;
    
    public MergeSplitConversion() { }

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
        
        gapRate = 1.0/gapSizeInput.get();
    }
    
    @Override
    public double proposal() {
        if (Randomizer.nextBoolean())
            return mergeProposal(); // Attempt merge
        else
            return splitProposal(); // Attempt split
    }
    
    double splitProposal() {
        double logHR = 0.0;
        
        if (acg.getConvCount()==0)
            return Double.NEGATIVE_INFINITY;
        
        // Select conversion
        Conversion conv = acg.getConversions()
                .get(Randomizer.nextInt(acg.getConvCount()));
        logHR -= Math.log(1.0/acg.getConvCount());
        
        if (conv.getSiteCount()<3)
            return Double.NEGATIVE_INFINITY;
        
        // Record original end of region:
        int origEnd = conv.getEndSite();
        
        // Select split point:
        int s = conv.getStartSite()
                + Randomizer.nextInt(conv.getSiteCount()-2);
        logHR -= Math.log(1.0/(double)(conv.getSiteCount()-2));
        
        // Select gap size
        int gap = (int)Randomizer.nextGeometric(gapRate);
        logHR -= gap*Math.log(1.0-gapRate) + Math.log(gapRate);
        
        if (s + gap + 2 > conv.getEndSite())
            return Double.NEGATIVE_INFINITY;
        
        // Select new departure height
        double depMin = conv.getNode1().getHeight();
        double depMax = conv.getNode1().getParent().getHeight();
        double depHeight = depMin + (depMax-depMin)*Randomizer.nextDouble();
        logHR -= Math.log(1.0/(depMax-depMin));
        
        // Select new arrival height
        double arrMin = conv.getNode2().getHeight();
        double arrMax;
        if (conv.getNode2().isRoot())
            arrMax = arrMin + 2.0*(conv.getHeight2() - arrMin);
        else
            arrMax = conv.getNode2().getParent().getHeight();            

        double arrHeight = arrMin + (arrMax-arrMin)*Randomizer.nextDouble();
        logHR -= Math.log(1.0/(arrMax-arrMin));
        
        if (arrHeight<depHeight)
            return Double.NEGATIVE_INFINITY;

        // Update original recombination
        conv.setEndSite(s);
        
        // Create new recombination
        Conversion newRecomb = new Conversion();
        newRecomb.setStartSite(s + gap + 2);
        newRecomb.setEndSite(origEnd);
        newRecomb.setNode1(conv.getNode1());
        newRecomb.setNode2(conv.getNode2());
        newRecomb.setHeight1(depHeight);
        newRecomb.setHeight2(arrHeight);
        acg.addConversion(newRecomb);
        
        // Include probability of reverse (merge) move in HR:
        logHR += Math.log(1.0/((double)(acg.getConvCount()-1)));
        
        return logHR;
    }
    
    double mergeProposal() {
        double logHR = 0.0;
        
        if (acg.getConvCount()<2)
            return Double.NEGATIVE_INFINITY;
        
        // Select a random pair of adjacent (on alignment) regions
        int r1idx = Randomizer.nextInt(acg.getConvCount()-1);
        int r2idx = r1idx + 1;
        Conversion recomb1 = acg.getConversions().get(r1idx);
        Conversion recomb2 = acg.getConversions().get(r2idx);
        logHR -= Math.log(1.0/((double)(acg.getConvCount()-1)));
        
        // Ensure a split operator could have generated this pair
        if (recomb1.getNode1() != recomb2.getNode1()
                || recomb1.getNode2() != recomb2.getNode2())
            return Double.NEGATIVE_INFINITY;
        
        // Record stuff needed for HR calculation:
        int gap = recomb2.getStartSite() - recomb1.getEndSite() - 2;
        double depRange = recomb1.getNode1().getParent().getHeight()
                - recomb1.getNode1().getHeight();
        double arrRange;
        if (recomb1.getNode2().isRoot())
            arrRange = 2.0*(recomb1.getHeight2()
                    - recomb1.getNode2().getHeight());
        else
            arrRange = recomb1.getNode2().getParent().getHeight()
                    - recomb1.getNode2().getHeight();
        
        // Final check to make sure height of recomb2's arrival point is
        // within the range available to the split move
        if (Math.abs(recomb2.getHeight2()-recomb1.getHeight2()) > arrRange)
            return Double.NEGATIVE_INFINITY;
        
        // Extend left-most region to cover entire region originally
        // spanned by both regions
        recomb1.setEndSite(recomb2.getEndSite());
        
        // Remove the right-most recombination
        acg.deleteConversion(recomb2);
        
        // Include probability of reverse (split) move in HR:
        logHR += Math.log(1.0/acg.getConvCount())
                + Math.log(1.0/((double)(recomb1.getSiteCount()-2)))
                + gap*Math.log(1.0-gapRate) + Math.log(gapRate)
                + Math.log(1.0/(depRange*arrRange));
        
        return logHR;
    }
}
