package argbeast.operators;


import argbeast.Recombination;
import argbeast.RecombinationGraph;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;

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

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which adds and removes recombinations to/from an ARG.")
public class AddRemoveRecombination extends RecombinationGraphOperator {
    
    public Input<RealParameter> rhoInput = new Input<RealParameter>("rho",
            "Recombination rate parameter.", Validate.REQUIRED);
    
    public Input<RealParameter> deltaInput = new Input<RealParameter>("delta",
            "Tract length parameter.", Validate.REQUIRED);
    
    private RecombinationGraph arg;
    
    public AddRemoveRecombination() { };
    
    @Override
    public void initAndValidate() {
        arg = argInput.get();
    };

    @Override
    public double proposal() {
        double logHGF = 0;

        if (Randomizer.nextDouble()<0.5) {
            
            // Add
            
            logHGF += Math.log(1.0/(arg.getNRecombs()+1));
            logHGF -= drawNewRecomb();
            
        } else {
            
            // Remove
            
            // Select recombination to remove:
            Recombination recomb = arg.getRecombinations().get(Randomizer.nextInt(arg.getNRecombs())+1);
            
            // Calculate HGF
            logHGF += getRecombProb(recomb);
            logHGF -= Math.log(1.0/arg.getNRecombs());
            
            // Remove recombination
            arg.deleteRecombination(recomb);
        }
        
        return logHGF;
    }
    
    private double drawNewRecomb() {
        double logP = 0;
        
        return logP;
    }
    
    private double getRecombProb(Recombination recomb) {
        double logP = 0;
        
        return logP;
    }
    
}
