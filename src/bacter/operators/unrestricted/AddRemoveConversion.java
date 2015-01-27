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
package bacter.operators.unrestricted;

import bacter.operators.EdgeCreationOperator;
import bacter.Conversion;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import feast.input.In;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which adds and removes conversions to/from an ACG.")
public class AddRemoveConversion extends EdgeCreationOperator {
    
    public Input<RealParameter> rhoInput = new In<RealParameter>("rho",
            "Conversion rate parameter.").setRequired();
    
    public Input<RealParameter> deltaInput = new In<RealParameter>("delta",
            "Tract length parameter.").setRequired();
    
    public AddRemoveConversion() { };
    
    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
    };

    @Override
    public double proposal() {
        double logHGF = 0;
        
        if (Randomizer.nextBoolean()) {
            
            // Add
            
            //logHGF += Math.log(1.0/(acg.getConvCount()+1));
            
            //logHGF -= drawNewConversion();
            drawNewConversion();
            
            if (!acg.isValid())
                return Double.NEGATIVE_INFINITY;
            
        } else {
            
            // Remove
            
            if (acg.getConvCount()==0)
                return Double.NEGATIVE_INFINITY;
            
            // Select conversion to remove:
            Conversion conv = acg.getConversions().get(
                    Randomizer.nextInt(acg.getConvCount()));
            
            // Calculate HGF
            //logHGF += getConversionProb(conv);
            //logHGF -= Math.log(1.0/acg.getConvCount());
            
            // Remove conversion
            acg.deleteConversion(conv);

        }

        return logHGF;
    }
    
    /**
     * Add new conversion to ACG, returning the probability density of the
     * new edge and converted region location.
     * 
     * @return log of proposal density 
     */
    public double drawNewConversion() {
        double logP = 0;

        Conversion newConversion = new Conversion();
        
        logP += attachEdge(newConversion);
        
        // Draw location of converted region.
        int startSite, endSite;
        double u = Randomizer.nextDouble()*(deltaInput.get().getValue() + acg.getSequenceLength());
        if (u<deltaInput.get().getValue()) {
            startSite = 0;
            logP += Math.log(deltaInput.get().getValue()
                /(deltaInput.get().getValue() + acg.getSequenceLength()));
        } else {
            startSite = (int)(u-deltaInput.get().getValue());
            logP += Math.log(1.0/(deltaInput.get().getValue()
                + acg.getSequenceLength()));
        }

        endSite = startSite + (int)Randomizer.nextGeometric(1.0/deltaInput.get().getValue());
        endSite = Math.min(endSite, acg.getSequenceLength()-1);

        // Probability of end site:
        double probEnd = Math.pow(1.0-1.0/deltaInput.get().getValue(),
            endSite-startSite)/ deltaInput.get().getValue();
        
        // Include probability of going past the end:
        if (endSite == acg.getSequenceLength()-1)
            probEnd += Math.pow(1.0-1.0/deltaInput.get().getValue(),
                    acg.getSequenceLength()-startSite);

        logP += Math.log(probEnd);

        newConversion.setStartSite(startSite);
        newConversion.setEndSite(endSite);

        acg.addConversion(newConversion);
        
        return logP;
    }
      
    /**
     * Obtain proposal density for the move which results in the current state
     * by adding the conversion conv to a state without that recombination.
     * 
     * @param conv
     * @return log of proposal density
     */
    public double getConversionProb(Conversion conv) {
        double logP = 0;
        
        logP += getEdgeAttachmentProb(conv);
        
        // Calculate probability of converted region.
        if (conv.getStartSite()==0)
            logP += Math.log((deltaInput.get().getValue() + 1)
                /(deltaInput.get().getValue() + acg.getSequenceLength()));
        else
            logP += Math.log(1.0/(deltaInput.get().getValue()
                + acg.getSequenceLength()));

        // Probability of end site:
        double probEnd = Math.pow(1.0-1.0/deltaInput.get().getValue(),
            conv.getEndSite() - conv.getStartSite())
            / deltaInput.get().getValue();
        
        // Include probability of going past the end:
        if (conv.getEndSite() == acg.getSequenceLength()-1)
            probEnd += Math.pow(1.0-1.0/deltaInput.get().getValue(),
                    acg.getSequenceLength()-conv.getStartSite());

        logP += Math.log(probEnd);
        
        return logP;
    }
}
