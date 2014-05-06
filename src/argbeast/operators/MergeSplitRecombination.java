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

package argbeast.operators;

import argbeast.Recombination;
import beast.core.Description;
import beast.core.Input;
import beast.util.Randomizer;
import feast.input.In;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which in one direction splits an existing recombination "
        + "in two and in the other merges two adjacent recombinations.")
public class MergeSplitRecombination extends RecombinationGraphOperator {
    
    public Input<Double> gapSizeInput = new In<Double>("expectedGapSize",
            "Expected gap size.").setDefault(10.0);

    private double gapRate;
    
    public MergeSplitRecombination() { }

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
        
        gapRate = 1.0/gapSizeInput.get();
    }
    
    

    @Override
    public double proposal() {
        
        if (Randomizer.nextBoolean()) {
            // Attempt merge
            
            return mergeProposal();
        } else {
            // Attempt split
            
            return splitProposal();
        }
    }
    
    private double splitProposal() {
        
        double logHR = 0.0;
        
        if (arg.getNRecombs()==0)
            return Double.NEGATIVE_INFINITY;
        
        // Select recombination
        Recombination recomb = arg.getRecombinations()
                .get(Randomizer.nextInt(arg.getNRecombs())+1);
        
        logHR -= Math.log(1.0/arg.getNRecombs());
        
        if (recomb.getSiteCount()<3)
            return Double.NEGATIVE_INFINITY;
        
        // Select split point:
        int s1 = recomb.getStartLocus()
                + Randomizer.nextInt(recomb.getSiteCount()-2) + 1;
        
        logHR -= Math.log(1.0/(double)(recomb.getSiteCount()-2));
        
        // Select gap size
        int gap = 1 + (int)Randomizer.nextGeometric(gapRate);
        
        logHR -= gap*Math.log(1-gapRate) + Math.log(gapRate);
        
        return 0;
    }
    
    private double mergeProposal() {
        return 0;
    }
}
