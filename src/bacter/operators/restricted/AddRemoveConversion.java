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
package bacter.operators.restricted;

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
    
    public Input<RealParameter> rhoInput = new Input<>("rho",
            "Conversion rate parameter.", Input.Validate.REQUIRED);
    
    public Input<RealParameter> deltaInput = new Input<>("delta",
            "Tract length parameter.", Input.Validate.REQUIRED);
    
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
            
            logHGF += Math.log(1.0/(acg.getTotalConvCount()+1));
            
            try {
                logHGF -= drawNewConversion();
            } catch (ProposalFailed ex) {
                return Double.NEGATIVE_INFINITY;
            }
            
            if (!acg.isValid())
                return Double.NEGATIVE_INFINITY;
            
        } else {
            
            // Remove
            
            if (acg.getTotalConvCount()==0)
                return Double.NEGATIVE_INFINITY;
            
            // Select conversion to remove:
            Conversion conv = chooseConversion();

            // Calculate HGF
            logHGF += getConversionProb(conv);
            logHGF -= Math.log(1.0/acg.getTotalConvCount());
            
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
    public double drawNewConversion() throws ProposalFailed {
        double logP = 0;

        Conversion newConversion = new Conversion();
        
        logP += attachEdge(newConversion);
        
        // Draw location of converted region.  Currently draws start locus 
        // uniformly from among available unconverted loci and draws the tract
        // length from an exponential distribution.  If the end of the tract
        // exceeds the start of the next region, the proposal is rejected.
        
        int convertableLength = getConvertableSiteCount(null);
        if (convertableLength==0)
            throw new ProposalFailed();
        
        int z = Randomizer.nextInt(convertableLength);
        logP += Math.log(1.0/convertableLength);
        
        for (Conversion recomb : acg.getConversions()) {
            
            if (z<recomb.getStartSite()-1)
                break;
            
            z += recomb.getEndSite()-Math.max(0,recomb.getStartSite()-1);
        }
        
        newConversion.setStartSite(z);
        
        int convertedLength = (int)Randomizer
                .nextGeometric(1.0/deltaInput.get().getValue());
        logP += convertedLength*Math.log(1.0-1.0/deltaInput.get().getValue())
                + Math.log(1.0/deltaInput.get().getValue());
                        
        newConversion.setEndSite(newConversion.getStartSite()+convertedLength);
        acg.addConversion(newConversion);

        // Abort if new conversion creates an invalid region map.
        if (!regionMapIsValid())
            throw new ProposalFailed();
        
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
        int convertableLength = getConvertableSiteCount(conv);
        if (convertableLength==0)
            return Double.NEGATIVE_INFINITY;
        
        logP += Math.log(1.0/convertableLength)
                + (conv.getEndSite()-conv.getStartSite())
                *Math.log(1.0-1.0/deltaInput.get().getValue())
                + Math.log(1.0/deltaInput.get().getValue());
        
        return logP;
    }
    
    /**
     * Get total number of sites available to convert, optionally skipping
     * one recombination.
     * 
     * @param skip recombination to skip (may be null)
     * @return convertible site count
     */
    private int getConvertableSiteCount(Conversion skip) {
        int count = 0;
        
        int l=0;
        for (int ridx=0; ridx<acg.getConvCount(); ridx++) {
            Conversion recomb = acg.getConversions().get(ridx);
            if (recomb == skip)
                continue;
            
            count += Math.max(0, recomb.getStartSite()-l+1);
            
            l = recomb.getEndSite()+2;
        }
        
        count += Math.max(0, acg.getSequenceLength() - l); // L - 1 - l + 1
        
        return count;
    }
}
