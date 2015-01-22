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

import bacter.operators.ACGOperator;
import bacter.Conversion;
import beast.core.Description;
import beast.core.Input;
import beast.util.Randomizer;
import feast.input.In;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which moves the alignment region affected "
        + "by a randomly-selected conversion event.")
public class ConvertedRegionShift extends ACGOperator {

    public Input<Double> apertureSizeInput = new In<Double>(
            "apertureSize",
            "Relative size (with respect to alignment size) of aperture "
                    + "within which new location of region edge is chosen "
                    + "uniformly. (Default 0.01, ie. 1%)").setDefault(0.01);

    public ConvertedRegionShift() { }

    @Override
    public double proposal() {
        
        if (acg.getConvCount()<1)
            return Double.NEGATIVE_INFINITY;
        
        int ridx = Randomizer.nextInt(acg.getConvCount());
        Conversion recomb = acg.getConversions().get(ridx);
        
        int radius = (int)Math.round(argInput.get().getSequenceLength()
                *apertureSizeInput.get())/2;

        int delta = Randomizer.nextInt(radius*2 + 1) - radius;
        
        if (delta>0) {
            int maxDelta;
            if (ridx<acg.getConvCount()-1)
                maxDelta = acg.getConversions().get(ridx+1).getStartSite() - 2
                        - recomb.getEndSite();
            else
                maxDelta = acg.getSequenceLength() - 1 - recomb.getEndSite();
            
            if (delta>maxDelta)
                return Double.NEGATIVE_INFINITY;
        } else {
            int minDelta;
            if (ridx>0)
                minDelta = acg.getConversions().get(ridx-1).getEndSite() + 2
                        - recomb.getStartSite();
            else
                minDelta = 0 - recomb.getStartSite();
            
            if (delta<minDelta)
                return Double.NEGATIVE_INFINITY;
        }
        
        recomb.setStartSite(recomb.getStartSite()+delta);
        recomb.setEndSite(recomb.getEndSite()+delta);
        
        return 0.0;
    }
    
}
