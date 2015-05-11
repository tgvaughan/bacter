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

import bacter.operators.EdgeCreationOperator;
import bacter.Conversion;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.util.Randomizer;
import feast.input.In;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Selects adjacent pair of conversions, deletes one, then expands "
        + "the other to encompass the region covered by the first.")
public class ReplaceConversion extends EdgeCreationOperator {
    
    public Input<RealParameter> rhoInput = new In<RealParameter>("rho",
            "Conversion rate parameter").setRequired();

    @Override
    public double proposal() {
        double logP = 0.0;

        double pRec = 1.0 - Math.exp(
                -0.5*acg.getClonalFrameLength()*rhoInput.get().getValue());

        Alignment alignment = chooseAlignment();
        
        if (Randomizer.nextBoolean()) {
            // Delete

            if (acg.getConvCount(alignment)<2)
                return Double.NEGATIVE_INFINITY;
        
            boolean right = Randomizer.nextBoolean();
            int ridx, gapSize;
            if (right) {
                ridx = Randomizer.nextInt(acg.getConvCount(alignment)-1);
                gapSize = acg.getConversions(alignment).get(ridx+1).getStartSite()
                        - acg.getConversions(alignment).get(ridx).getEndSite() - 1;
            } else {
                ridx = Randomizer.nextInt(acg.getConvCount(alignment)-1) + 1;
                gapSize = acg.getConversions(alignment).get(ridx).getStartSite()
                        - acg.getConversions(alignment).get(ridx-1).getEndSite() - 1;
            }
            
            Conversion recomb = acg.getConversions(alignment).get(ridx);
            
            logP -= Math.log(1.0/(acg.getConvCount(alignment)-1));
            
            if (right) {
                logP += getEdgeAttachmentProb(acg.getConversions(alignment).get(ridx+1));
                
                recomb.setEndSite(acg.getConversions(alignment).get(ridx+1).getEndSite());
                acg.deleteConversion(acg.getConversions(alignment).get(ridx+1));

            } else {
                logP += getEdgeAttachmentProb(acg.getConversions(alignment).get(ridx-1));
                
                recomb.setStartSite(acg.getConversions(alignment).get(ridx-1).getStartSite());
                acg.deleteConversion(acg.getConversions(alignment).get(ridx-1));
            }
            
            logP += Math.log(1.0/acg.getConvCount(alignment))
                    + Math.log(1.0/(recomb.getEndSite()-recomb.getStartSite()-1))
                    + (gapSize-1)*Math.log(1.0-pRec) + Math.log(pRec);
            
        } else {
            // Add
            
            if (acg.getConvCount(alignment)<1)
                return Double.NEGATIVE_INFINITY;
            
            boolean right = Randomizer.nextBoolean();
            int ridx = Randomizer.nextInt(acg.getConvCount(alignment));
            
            logP -= Math.log(1.0/acg.getConvCount(alignment));
            
            Conversion conv = acg.getConversions(alignment).get(ridx);
            
            if (conv.getEndSite()-conv.getStartSite() + 1 < 3)
                return Double.NEGATIVE_INFINITY;

            int gapSize = 1 + (int)Randomizer.nextGeometric(pRec);
            
            int startGap = conv.getStartSite() + 1 + Randomizer.nextInt(
                    conv.getEndSite() - conv.getStartSite() - 1);
            int endGap = startGap + gapSize - 1;
            
            if (endGap>=conv.getEndSite())
                return Double.NEGATIVE_INFINITY;
            
            logP -= (gapSize-1)*Math.log(1.0-pRec) + Math.log(pRec)
                    + Math.log(1.0/(conv.getEndSite()-conv.getStartSite()-1));
            
            Conversion newConv = new Conversion(conv.getAlignment());
            if (right) {
                newConv.setStartSite(endGap+1);
                newConv.setEndSite(conv.getEndSite());
                conv.setEndSite(startGap-1);
            } else {
                newConv.setStartSite(conv.getStartSite());
                newConv.setEndSite(startGap-1);
                conv.setStartSite(endGap+1);
            }
            logP -= attachEdge(newConv);
            
            acg.addConversion(newConv);
            
            logP += Math.log(1.0/(acg.getConvCount(alignment)-1));
        }
        
        return logP;
    }

}
