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

        // Total effective number of possible start sites
        double alpha = acg.getTotalSequenceLength()
                + acg.getAlignments().size()*deltaInput.get().getValue();

        // Draw location of converted region.
        int startSite = -1;
        int endSite = -1;
        Alignment alignment = null;

        double u = Randomizer.nextDouble()*alpha;
        for (Alignment thisAlignment : acg.getAlignments()) {
            if (u < deltaInput.get().getValue() + acg.getSequenceLength(thisAlignment)) {
                alignment = thisAlignment;

                if (u < deltaInput.get().getValue()) {
                    startSite = 0;
                    logP += Math.log(deltaInput.get().getValue() / alpha);
                } else {
                    startSite = (int) (u - deltaInput.get().getValue());
                    logP += Math.log(1.0 / alpha);
                }

                break;
            }

            u -= deltaInput.get().getValue() + acg.getSequenceLength(thisAlignment);
        }

        if (alignment == null)
            throw new IllegalStateException("Programmer error: " +
                    "loop in drawAffectedRegion() fell through.");

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

        conv.setAlignment(alignment);
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

        // Total effective number of possible start sites
        double alpha = acg.getTotalSequenceLength()
                + acg.getAlignments().size()*deltaInput.get().getValue();


        // Calculate probability of converted region.
        if (conv.getStartSite()==0)
            logP += Math.log((deltaInput.get().getValue() + 1) / alpha);
        else
            logP += Math.log(1.0 / alpha);

        // Probability of end site:
        double probEnd = Math.pow(1.0-1.0/deltaInput.get().getValue(),
            conv.getEndSite() - conv.getStartSite())
            / deltaInput.get().getValue();
        
        // Include probability of going past the end:
        if (conv.getEndSite() == acg.getSequenceLength(conv.getAlignment())-1)
            probEnd += Math.pow(1.0-1.0/deltaInput.get().getValue(),
                    acg.getSequenceLength(conv.getAlignment())-conv.getStartSite());

        logP += Math.log(probEnd);

        return logP;
    }
}
