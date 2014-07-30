
import bacter.Recombination;
import bacter.operators.RecombinationGraphOperator;
import beast.core.Description;
import beast.core.Input;
import beast.util.Randomizer;
import feast.input.In;

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

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which selects a pair of adjacent regions and jointly "
        + "shifts the closest edge pair, keeping the distance between them "
        + "the same.")
public class ConvertedRegionPairBoundaryShift extends RecombinationGraphOperator{

    public Input<Double> apertureSizeInput = new In<Double>(
            "apertureSize",
            "Relative size (with respect to alignment size) of aperture "
            + "within which new location of region edge is chosen "
            + "uniformly. (Default 0.01, ie. 1%)").setDefault(0.01);
    
    public ConvertedRegionPairBoundaryShift() { }
    
    @Override
    public double proposal() {
        if (arg.getNRecombs()<2)
            return Double.NEGATIVE_INFINITY;
        
        // Select random recombination to be left of pair:
        int ridx = Randomizer.nextInt(arg.getNRecombs()-1) + 1;
        Recombination leftRecomb = arg.getRecombinations().get(ridx);
        Recombination rightRecomb = arg.getRecombinations().get(ridx+1);
        
        int radius = (int)Math.round(arg.getSequenceLength()
                *apertureSizeInput.get())/2;
        
        int minOffset = -(leftRecomb.getEndSite() - leftRecomb.getStartSite());
        int maxOffset = rightRecomb.getEndSite() - rightRecomb.getStartSite();
        
        int offset = Randomizer.nextInt(2*radius + 1) - radius;
        
        if (offset<minOffset || offset>maxOffset)
            return Double.NEGATIVE_INFINITY;
        
        leftRecomb.setEndSite(leftRecomb.getEndSite()+offset);
        rightRecomb.setStartSite(rightRecomb.getStartSite()+offset);
        
        return 0;
    }
    
}
