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

package bacter.operators;

import bacter.Conversion;
import beast.core.Description;
import beast.core.Input;
import beast.util.Randomizer;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which moves one edge of the alignment region affected "
        + "by a randomly-selected conversion event.")
public class ConvertedRegionBoundaryShift extends ACGOperator {

    public Input<Double> apertureSizeInput = new Input<>(
            "apertureSize",
            "Relative size (with respect to alignment size) of aperture "
                    + "within which new location of region edge is chosen "
                    + "uniformly. (Default 0.01, ie. 1%)", 0.01);

    public ConvertedRegionBoundaryShift() { }

    @Override
    public double proposal() {
        
        if (acg.getTotalConvCount()<1 || acg.wholeLocusModeOn())
            return Double.NEGATIVE_INFINITY;
        
        // Select random conversion and region edge:
        Conversion conv = chooseConversion();
        boolean moveStart = Randomizer.nextBoolean();
        
        int currentLocus, minLocus, maxLocus;
        if (moveStart) {
            currentLocus = conv.getStartSite();
            maxLocus = conv.getEndSite();
            minLocus = 0;
        } else {
            currentLocus = conv.getEndSite();
            minLocus = conv.getStartSite();
            maxLocus = conv.getLocus().getSiteCount()-1;
        }
        
        int radius = (int)Math.round(conv.getLocus().getSiteCount()
                * apertureSizeInput.get())/2;
        
        int newLocus = currentLocus + Randomizer.nextInt(2*radius+1)-radius;
        
        if (newLocus < minLocus || newLocus > maxLocus)
            return Double.NEGATIVE_INFINITY;

        if (moveStart)
            conv.setStartSite(newLocus);
        else
            conv.setEndSite(newLocus);

        assert !acg.isInvalid() : "CRBS produced invalid state.";
        
        return 0;
    }
    
}
