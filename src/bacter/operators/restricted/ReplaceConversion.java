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
        
        if (Randomizer.nextBoolean()) {
            // Delete

            if (acg.getConvCount()<2)
                return Double.NEGATIVE_INFINITY;
        
            boolean right = Randomizer.nextBoolean();
            int ridx, gapSize;
            if (right) {
                ridx = Randomizer.nextInt(acg.getConvCount()-1);
                gapSize = acg.getConversions().get(ridx+1).getStartSite()
                        - acg.getConversions().get(ridx).getEndSite() - 1;
            } else {
                ridx = Randomizer.nextInt(acg.getConvCount()-1) + 1;
                gapSize = acg.getConversions().get(ridx).getStartSite()
                        - acg.getConversions().get(ridx-1).getEndSite() - 1;
            }
            
            Conversion recomb = acg.getConversions().get(ridx);
            
            logP -= Math.log(1.0/(acg.getConvCount()-1));
            
            if (right) {
                logP += getEdgeAttachmentProb(acg.getConversions().get(ridx+1));
                
                recomb.setEndSite(acg.getConversions().get(ridx+1).getEndSite());
                acg.deleteConversion(acg.getConversions().get(ridx+1));

            } else {
                logP += getEdgeAttachmentProb(acg.getConversions().get(ridx-1));
                
                recomb.setStartSite(acg.getConversions().get(ridx-1).getStartSite());
                acg.deleteConversion(acg.getConversions().get(ridx-1));
            }
            
            logP += Math.log(1.0/acg.getConvCount())
                    + Math.log(1.0/(recomb.getEndSite()-recomb.getStartSite()-1))
                    + (gapSize-1)*Math.log(1.0-pRec) + Math.log(pRec);
            
        } else {
            // Add
            
            if (acg.getConvCount()<1)
                return Double.NEGATIVE_INFINITY;
            
            boolean right = Randomizer.nextBoolean();
            int ridx = Randomizer.nextInt(acg.getConvCount());
            
            logP -= Math.log(1.0/acg.getConvCount());
            
            Conversion recomb = acg.getConversions().get(ridx);
            
            if (recomb.getEndSite()-recomb.getStartSite() + 1 < 3)
                return Double.NEGATIVE_INFINITY;

            int gapSize = 1 + (int)Randomizer.nextGeometric(pRec);
            
            int startGap = recomb.getStartSite() + 1 + Randomizer.nextInt(
                    recomb.getEndSite() - recomb.getStartSite() - 1);
            int endGap = startGap + gapSize - 1;
            
            if (endGap>=recomb.getEndSite())
                return Double.NEGATIVE_INFINITY;
            
            logP -= (gapSize-1)*Math.log(1.0-pRec) + Math.log(pRec)
                    + Math.log(1.0/(recomb.getEndSite()-recomb.getStartSite()-1));
            
            Conversion newConv = new Conversion();
            if (right) {
                newConv.setStartSite(endGap+1);
                newConv.setEndSite(recomb.getEndSite());
                recomb.setEndSite(startGap-1);
            } else {
                newConv.setStartSite(recomb.getStartSite());
                newConv.setEndSite(startGap-1);
                recomb.setStartSite(endGap+1);
            }
            logP -= attachEdge(newConv);
            
            acg.addConversion(newConv);
            
            logP += Math.log(1.0/(acg.getConvCount()-1));
        }
        
        return logP;
    }

}
