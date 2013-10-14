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
package argbeast;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.Coalescent;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Appoximation to the coalescent with gene conversion.")
public class GCCoalescentApprox extends Coalescent {
    
    public Input<RealParameter> rhoInput = new Input<RealParameter>("rho",
            "Recombination rate parameter.", Validate.REQUIRED);
    
    public Input<RealParameter> deltaInput = new Input<RealParameter>("delta",
            "Tract length parameter.", Validate.REQUIRED);
    
    RecombinationGraph arg;
    int sequenceLength;
    
    @Override
    public void initAndValidate() throws Exception {
        if (!(treeInput.get() instanceof RecombinationGraph))
            throw new Exception("Error: GCCoalescentApprox only applicable to"
                    + " RecombinationGraph objects.");
        
        arg = (RecombinationGraph)treeInput.get();
        sequenceLength = arg.getSequenceLength();
        
        super.initAndValidate();
    }
    
    @Override
    public double calculateLogP() throws Exception {
        
        logP = clonalFrameLogP()
                + conditionalCoalescentLogP()
                + convertedRegionLogP();
        
        return logP;        
    }
    
    /**
     * Calculate probability of clonal frame under the standard coalescent.
     * 
     * @return log(P)
     */
    double clonalFrameLogP() {
        return calculateLogLikelihood(treeIntervalsInput.get(),
                popSizeInput.get());
    }
    
    /**
     * Compute probability of recombinant edges under conditional coalescent.
     * @return log(P)
     */
    double conditionalCoalescentLogP() {
        return 0.0;
    }
    
    /**
     * Compute probability of number and genome extent of converted segments.
     * @return log(P)
     */
    double convertedRegionLogP() {
        
        double crLogP = 0.0;
        
        List<Recombination> recombs = arg.getRecombinations();
        
        double rho = rhoInput.get().getValue();
        double deltaInv = 1.0/deltaInput.get().getValue();
        double cfLength = arg.getClonalFrameLength();
        
        double alpha = 0.5*rho*cfLength/sequenceLength;
        
        // Account for starting state
        double pStartCF = 1.0/(alpha/deltaInv+1);        
        if (recombs.size()>0 && recombs.get(0).startLocus>0)
            crLogP += Math.log(pStartCF);
        else
            crLogP += Math.log(1.0 - pStartCF);
        
        for (Recombination recomb : recombs) {
            if (recomb.startLocus>0)
                crLogP += Math.log(alpha);
            
            crLogP += (recomb.endLocus-recomb.startLocus)*Math.log(1-deltaInv);
            
            if (recomb.endLocus<sequenceLength-1)
                crLogP += Math.log(deltaInv);
        }
        
        return crLogP;
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }
}
