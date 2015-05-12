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
        
        if (acg.getTotalConvCount()<1)
            return Double.NEGATIVE_INFINITY;
        
        Conversion conv = chooseConversion();
        
        int radius = (int)Math.round(argInput.get().getSequenceLength(conv.getAlignment())
            *apertureSizeInput.get())/2;

        int delta = Randomizer.nextInt(radius*2 + 1) - radius;
        
        if (conv.getEndSite() + delta > acg.getSequenceLength(conv.getAlignment()) - 1)
            return Double.NEGATIVE_INFINITY;

        if (conv.getStartSite() + delta<0)
            return Double.NEGATIVE_INFINITY;
        
        conv.setStartSite(conv.getStartSite()+delta);
        conv.setEndSite(conv.getEndSite()+delta);
        
        return 0.0;
    }
    
}
