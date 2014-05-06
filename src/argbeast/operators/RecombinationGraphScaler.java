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

package argbeast.operators;

import argbeast.Recombination;
import argbeast.RecombinationGraph;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import feast.input.In;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Scaling operator for recombination graphs.")
public class RecombinationGraphScaler extends RecombinationGraphOperator {
    
    public Input<List<RealParameter>> parametersInput = new In<List<RealParameter>>(
            "parameter", "Parameter to scale with ARG.")
            .setDefault(new ArrayList<RealParameter>());

    public Input<List<RealParameter>> parametersInverseInput = new In<List<RealParameter>>(
            "parameterInverse", "Parameter to scale inversely with ARG.")
            .setDefault(new ArrayList<RealParameter>());
    
    public Input<Double> scaleParamInput = new In<Double>("scaleFactor",
            "Scale factor tuning parameter.  Must be < 1.").setRequired();
    
    private double scaleParam;
    
    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
        scaleParam = scaleParamInput.get();
    }
    
    @Override
    public double proposal() {
        
        // Keep track of number of positively scaled elements minus
        // negatively scaled elements.
        int count = 0;

        // Choose scaling factor:
        double f = scaleParam + Randomizer.nextDouble()*(1.0/scaleParam - scaleParam);
        
        // Scale clonal frame:
        for (Node node : arg.getInternalNodes()) {
            node.setHeight(node.getHeight()*f);
            count += 1;
        }
        
        // Scale recombinant edges:
        for (Recombination recomb : arg.getRecombinations()) {
            if (recomb == null)
                continue;
            
            recomb.setHeight1(recomb.getHeight1()*f);
            recomb.setHeight2(recomb.getHeight2()*f);
            count += 2;
            
            if (recomb.getHeight1()<recomb.getNode1().getHeight()
                    || recomb.getHeight2()<recomb.getNode2().getHeight())
                return Double.NEGATIVE_INFINITY;
        }
        
        // Check for illegal node heights:
        for (Node node : arg.getExternalNodes()) {
            if (node.getHeight()>node.getParent().getHeight())
                return Double.NEGATIVE_INFINITY;
        }
        
        // Scale parameters
        
        for (RealParameter param : parametersInput.get()) {
            try {
                param.scale(f);
            } catch (Exception e) {
                
                // Scale throws Exception if param has been scaled outside its
                // bounds.  Needs to change!
                return Double.NEGATIVE_INFINITY;
            }
            
            count += param.getDimension();
        }
        
        for (RealParameter paramInv : parametersInverseInput.get()) {
            try {
                paramInv.scale(1.0/f);
            } catch (Exception e) {
                
                // Scale throws Exception if param has been scaled outside its
                // bounds.  Needs to change!
                return Double.NEGATIVE_INFINITY;
            }
            
            count -= paramInv.getDimension();
        }
        
        // Return log of Hastings ratio:
        return (count-2)*Math.log(f);
    }
    
}
