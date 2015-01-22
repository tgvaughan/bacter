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

package bacter.operators;

import bacter.operators.ConversionGraphOperator;
import bacter.Conversion;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import feast.input.In;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which moves connection points of recombinant edges "
        + "about on clonal frame.")
public class ConvertedEdgeSlide extends ConversionGraphOperator {

    public Input<Double> scaleBoundInput = new In<Double>("scaleBound",
            "Determines bounds of height scaling: [1/scaleBound, scaleBound]. "
                    + "Default is 0.8.").setDefault(0.8);

    private double scaleMin, scaleMax;

    public ConvertedEdgeSlide() { }

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
        scaleMin = Math.min(scaleBoundInput.get(), 1.0/scaleBoundInput.get());
        scaleMax = 1.0/scaleMin;
    }
    
    @Override
    public double proposal() {
        
        double logHR = 0.0;
        
        if (acg.getConvCount()==0)
            return Double.NEGATIVE_INFINITY;

        // Select edge at random:
        Conversion recomb = acg.getConversions().get(
                Randomizer.nextInt(acg.getConvCount()));
        
        // Decide whether to move departure or arrival point
        boolean moveDeparture = Randomizer.nextBoolean();
        
        // Get current (old) attachment height
        double oldHeight;
        if (moveDeparture) {
            oldHeight = recomb.getHeight1();
        } else {
            oldHeight = recomb.getHeight2();
        }
        
        // Choose scale factor
        double f = Randomizer.nextDouble()*(scaleMax - scaleMin) + scaleMin;
        
        // Set new height
        double newHeight = f*oldHeight;
        
        // Check for boundary violation
        if (moveDeparture) {
            if (newHeight>recomb.getHeight2() || newHeight>acg.getRoot().getHeight())
                return Double.NEGATIVE_INFINITY;
        } else {
            if (newHeight<recomb.getHeight1())
                return Double.NEGATIVE_INFINITY;
        }
        
        // Scale factor HR contribution
        logHR += -Math.log(f);

        // Get node below current (old) attachment point
        Node nodeBelow;
        if (moveDeparture)
            nodeBelow = recomb.getNode1();
        else
            nodeBelow = recomb.getNode2();

        // Choose node below new attachment point
        if (f<1.0) {
            while (newHeight<nodeBelow.getHeight()) {
                if (nodeBelow.isLeaf())
                    return Double.NEGATIVE_INFINITY;
                
                if (Randomizer.nextBoolean())
                    nodeBelow = nodeBelow.getLeft();
                else
                    nodeBelow = nodeBelow.getRight();
                
                logHR += -Math.log(0.5);
            }
        } else {
            while (!nodeBelow.isRoot() && newHeight>nodeBelow.getParent().getHeight()) {
                nodeBelow = nodeBelow.getParent();
                logHR += Math.log(0.5);
            }
        }
        
        // Write changes back to recombination object
        if (moveDeparture) {
            recomb.setHeight1(newHeight);
            recomb.setNode1(nodeBelow);
        } else {
            recomb.setHeight2(newHeight);
            recomb.setNode2(nodeBelow);
        }
        
        return logHR;
    }
    
}
