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

import bacter.Recombination;
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
public class ReplaceRecombination extends EdgeCreationOperator {
    
    public Input<RealParameter> rhoInput = new In<RealParameter>("rho",
            "Conversion rate parameter").setRequired();

    @Override
    public double proposal() {
        double logP = 0.0;
        
        double pRec = 1.0 - Math.exp(
                -0.5*arg.getClonalFrameLength()*rhoInput.get().getValue());
        
        if (Randomizer.nextBoolean()) {
            // Delete

            if (arg.getNRecombs()<2)
                return Double.NEGATIVE_INFINITY;
        
            boolean right = Randomizer.nextBoolean();
            int ridx, gapSize;
            if (right) {
                ridx = Randomizer.nextInt(arg.getNRecombs()-1) + 1;
                gapSize = arg.getRecombinations().get(ridx+1).getStartSite()
                        - arg.getRecombinations().get(ridx).getEndSite() - 1;
            } else {
                ridx = Randomizer.nextInt(arg.getNRecombs()-1) + 2;
                gapSize = arg.getRecombinations().get(ridx).getStartSite()
                        - arg.getRecombinations().get(ridx-1).getEndSite() - 1;
            }
            
            Recombination recomb = arg.getRecombinations().get(ridx);
            
            logP -= Math.log(1.0/(arg.getNRecombs()-1));
            
            if (right) {
                logP += getEdgeAttachmentProb(arg.getRecombinations().get(ridx+1));
                
                recomb.setEndSite(arg.getRecombinations().get(ridx+1).getEndSite());
                arg.deleteRecombination(arg.getRecombinations().get(ridx+1));

            } else {
                logP += getEdgeAttachmentProb(arg.getRecombinations().get(ridx-1));
                
                recomb.setStartSite(arg.getRecombinations().get(ridx-1).getStartSite());
                arg.deleteRecombination(arg.getRecombinations().get(ridx-1));
            }
            
            logP += Math.log(1.0/arg.getNRecombs())
                    + Math.log(1.0/(recomb.getEndSite()-recomb.getStartSite()-1))
                    + (gapSize-1)*Math.log(1.0-pRec) + Math.log(pRec);
            
        } else {
            // Add
            
            if (arg.getNRecombs()<1)
                return Double.NEGATIVE_INFINITY;
            
            boolean right = Randomizer.nextBoolean();
            int ridx = Randomizer.nextInt(arg.getNRecombs()) + 1;
            
            logP -= Math.log(1.0/arg.getNRecombs());
            
            Recombination recomb = arg.getRecombinations().get(ridx);
            
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
            
            Recombination newRecomb = new Recombination();
            if (right) {
                newRecomb.setStartSite(endGap+1);
                newRecomb.setEndSite(recomb.getEndSite());
                recomb.setEndSite(startGap-1);
            } else {
                newRecomb.setStartSite(recomb.getStartSite());
                newRecomb.setEndSite(startGap-1);
                recomb.setStartSite(endGap+1);
            }
            logP -= attachEdge(newRecomb);
            
            arg.addRecombination(newRecomb);
            
            logP += Math.log(1.0/(arg.getNRecombs()-1));
        }
        
        return logP;
    }

}
