/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
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
import beast.util.Randomizer;

/**
 *
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Performs a combined add/remove conversion operation, keeping " +
        "the total number of conversions constant.")
public class ReplaceConversion extends ConversionCreationOperator {

    @Override
    public double proposal() {
        if (acg.getTotalConvCount()==0)
            return Double.NEGATIVE_INFINITY;

        double logHGF = 0;

        // Select conversion to replace:
        Conversion conv = chooseConversion();

        // Incorporate conversion prior probability into HGF
        logHGF += getEdgeAttachmentProb(conv) + getAffectedRegionProb(conv);

        // Remove conversion
        acg.deleteConversion(conv);

        // Draw replacement conversion from prior, incoroporating
        // probability into HGF
        conv = new Conversion();
        logHGF -= attachEdge(conv) + drawAffectedRegion(conv);

        // Add replacement conversion
        acg.addConversion(conv);

        assert !acg.isInvalid() : "ReplaceConversion produced invalid state.";

        return logHGF;
    }
}
