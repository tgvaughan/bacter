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
import beast.evolution.tree.coalescent.Coalescent;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Appoximation to the coalescent with gene conversion.")
public class GCCoalescentApprox extends Coalescent {
    
    RecombinationGraph arg;
    
    @Override
    public void initAndValidate() throws Exception {
        if (!(treeInput.get() instanceof RecombinationGraph))
            throw new Exception("Error: GCCoalescentApprox only applicable to"
                    + " RecombinationGraph objects.");
        
        arg = (RecombinationGraph)treeInput.get();
        
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
        return 0.0;
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }
}
