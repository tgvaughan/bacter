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
import argbeast.RecombinationGraph;
import beast.core.Description;
import beast.core.Input;
import beast.util.Randomizer;
import feast.input.In;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which moves one edge of the alignment region affected "
        + "by a randomly-selected conversion event.")
public class ConvertedRegionBoundaryShift extends RecombinationGraphOperator {

    public Input<Double> apertureSizeInput = new In<Double>(
            "apertureSize",
            "Relative size (with respect to alignment size) of aperture "
                    + "within which new location of region edge is chosen "
                    + "uniformly. (Default 0.01, ie. 1%)").setDefault(0.01);

    public ConvertedRegionBoundaryShift() { }

    @Override
    public double proposal() {
        
        if (arg.getNRecombs()<1)
            return Double.NEGATIVE_INFINITY;
        
        // Select random recombination and region edge:
        int z = Randomizer.nextInt(arg.getNRecombs()*2);
        int ridx = z/2 + 1;
        Recombination recomb = arg.getRecombinations().get(ridx);
        boolean moveStart = (z%2 == 0);
        
        int currentLocus, minLocus, maxLocus;
        if (moveStart) {
            currentLocus = recomb.getStartSite();
            maxLocus = recomb.getEndSite();
            
            if (ridx>1)
                minLocus = arg.getRecombinations().get(ridx-1).getEndSite() + 2;
            else
                minLocus = 0;
        } else {
            currentLocus = recomb.getEndSite();
            minLocus = recomb.getStartSite();
            
            if (ridx<arg.getNRecombs())
                maxLocus = arg.getRecombinations().get(ridx+1).getStartSite() - 2;
            else
                maxLocus = arg.getSequenceLength()-1;
        }
        
        int radius = (int)Math.round(argInput.get().getSequenceLength()
                *apertureSizeInput.get())/2;
        
        int newLocus = currentLocus + Randomizer.nextInt(2*radius+1)-radius;
        
        if (newLocus < minLocus || newLocus > maxLocus)
            return Double.NEGATIVE_INFINITY;

        if (moveStart)
            recomb.setStartSite(newLocus);
        else
            recomb.setEndSite(newLocus);
        
        return 0;
    }
    
}