/*
 * Copyright (C) 2013 Tim Vaughan <tgvaughan@gmail.com>
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

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which moves one edge of the alignment region affected "
        + "by a randomly-selected conversion event.")
public class ConvertedRegionEdgeShift extends RecombinationGraphOperator {

    public Input<Double> apertureSizeInput = new Input<Double>(
            "apertureSize",
            "Relative size (with respect to alignment size) of aperture "
                    + "within which new location of region edge is chosen "
                    + "uniformly. (Default 0.01, ie. 1%)", 0.01);

    public ConvertedRegionEdgeShift() { }

    @Override
    public void initAndValidate() throws Exception { }
    
    @Override
    public double proposal() {

        if (argInput.get().getNRecombs()<1)
            return Double.NEGATIVE_INFINITY;
        
        // Select random recombination and region edge:
        int z = Randomizer.nextInt(argInput.get().getNRecombs()*2);
        int ridx = z/2 + 1;
        Recombination recomb = argInput.get().getRecombinations().get(ridx);
        boolean moveStart = (z%2 == 0);
        
        long currentLocus;
        if (moveStart)
            currentLocus = recomb.getStartLocus();
        else
            currentLocus = recomb.getEndLocus();
        
        // Identify boundaries of aperture
        long delta = Math.round(argInput.get().getSequenceLength()
                *apertureSizeInput.get())/2;
        
        long lower, upper;
        if (moveStart) {
            if (ridx==1)
                lower = 0;
            else
                lower = argInput.get().getRecombinations().get(ridx-1).getEndLocus()+2;
        } else
            lower = recomb.getStartLocus();
        
        lower = Math.max(lower, currentLocus-delta);
        
        if (moveStart) {
            upper = recomb.getEndLocus();
        } else {
            if (ridx==argInput.get().getNRecombs())
                upper = argInput.get().getSequenceLength()-1;
            else
                upper = argInput.get().getRecombinations().get(ridx+1).getStartLocus()-2;
        }
        upper = Math.min(upper, currentLocus+delta);
        
        // Select new site
        long locus = lower + Randomizer.nextInt((int)(upper-lower+1));
                
        // Perform shift
        if (moveStart)
            recomb.setStartLocus(locus);
        else
            recomb.setEndLocus(locus);
        
        // Calculate new boundaries for HR calculation
        long lowerPrime, upperPrime;
        if (moveStart) {
            if (ridx==1)
                lowerPrime = 0;
            else
                lowerPrime = argInput.get().getRecombinations().get(ridx-1).getEndLocus()+2;
        } else
            lowerPrime = recomb.getStartLocus();
        
        lowerPrime = Math.max(lowerPrime, currentLocus-delta);
        
        if (moveStart) {
            upperPrime = recomb.getEndLocus();
        } else {
            if (ridx==argInput.get().getNRecombs())
                upperPrime = argInput.get().getSequenceLength()-1;
            else
                upperPrime = argInput.get().getRecombinations().get(ridx+1).getStartLocus()-2;
        }
        upperPrime = Math.min(upperPrime, currentLocus+delta);

        return Math.log(((double)(upper-lower+1))/((double)(upperPrime-lowerPrime+1)));
    }
    
}
