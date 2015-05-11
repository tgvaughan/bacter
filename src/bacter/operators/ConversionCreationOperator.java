/*
 * Copyright (C) 2015 Tim Vaughan (tgvaughan@gmail.com)
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
import bacter.operators.EdgeCreationOperator;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.util.Randomizer;
import feast.input.In;

/**
 * Abstract class of ACG operators that use the clonal origin model as the
 * basis for adding new converted edges and their affected sites to an
 * existing ConversionGraph.
 *
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
public abstract class ConversionCreationOperator extends EdgeCreationOperator {

    public Input<RealParameter> deltaInput = new In<RealParameter>("delta",
            "Tract length parameter.").setRequired();
     
    /**
     * Choose region to be affected by this conversion.
     * 
     * @param conv Conversion object where these sites are stored.
     * @return log probability density of chosen attachment.
     */
    public double drawAffectedRegion(Conversion conv) {
        double logP = 0.0;

        Alignment alignment = conv.getAlignment();

        // Draw location of converted region.
        int startSite, endSite;
        double u = Randomizer.nextDouble()*(deltaInput.get().getValue() + acg.getSequenceLength(alignment));
        if (u<deltaInput.get().getValue()) {
            startSite = 0;
            logP += Math.log(deltaInput.get().getValue()
                /(deltaInput.get().getValue() + acg.getSequenceLength(alignment)));
        } else {
            startSite = (int)(u-deltaInput.get().getValue());
            logP += Math.log(1.0/(deltaInput.get().getValue()
                + acg.getSequenceLength(alignment)));
        }

        endSite = startSite + (int)Randomizer.nextGeometric(1.0/deltaInput.get().getValue());
        endSite = Math.min(endSite, acg.getSequenceLength(alignment)-1);

        // Probability of end site:
        double probEnd = Math.pow(1.0-1.0/deltaInput.get().getValue(),
            endSite-startSite)/ deltaInput.get().getValue();
        
        // Include probability of going past the end:
        if (endSite == acg.getSequenceLength(alignment)-1)
            probEnd += Math.pow(1.0-1.0/deltaInput.get().getValue(),
                    acg.getSequenceLength(alignment)-startSite);

        logP += Math.log(probEnd);

        conv.setStartSite(startSite);
        conv.setEndSite(endSite);

        return logP;
    }
    
    /**
     * Calculate probability of choosing region affected by the given
     * conversion under the ClonalOrigin model.
     * 
     * @param conv
     * @return log probability density
     */
    public double getAffectedRegionProb(Conversion conv) {
        double logP = 0.0;

        Alignment alignment = conv.getAlignment();

        // Calculate probability of converted region.
        if (conv.getStartSite()==0)
            logP += Math.log((deltaInput.get().getValue() + 1)
                /(deltaInput.get().getValue() + acg.getSequenceLength(alignment)));
        else
            logP += Math.log(1.0/(deltaInput.get().getValue()
                + acg.getSequenceLength(alignment)));

        // Probability of end site:
        double probEnd = Math.pow(1.0-1.0/deltaInput.get().getValue(),
            conv.getEndSite() - conv.getStartSite())
            / deltaInput.get().getValue();
        
        // Include probability of going past the end:
        if (conv.getEndSite() == acg.getSequenceLength(alignment)-1)
            probEnd += Math.pow(1.0-1.0/deltaInput.get().getValue(),
                    acg.getSequenceLength(alignment)-conv.getStartSite());

        logP += Math.log(probEnd);

        return logP;
    }
}
