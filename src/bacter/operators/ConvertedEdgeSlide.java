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

import bacter.Conversion;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which moves connection points of recombinant edges "
        + "about on clonal frame.")
public class ConvertedEdgeSlide extends ACGOperator {

    public Input<Double> apertureSizeInput = new Input<>("apertureSize",
            "Window size as a fraction of the clonal frame tree height."
                    + "Default is 0.1.", 0.1);

    public ConvertedEdgeSlide() { }

    @Override
    public void initAndValidate() {
        super.initAndValidate();
    }
    
    @Override
    public double proposal() {
        
        double logHR = 0.0;
        
        if (acg.getTotalConvCount()==0)
            return Double.NEGATIVE_INFINITY;

        // Select edge at random:
        Conversion conv = chooseConversion();
        
        // Decide whether to move departure or arrival point
        boolean moveDeparture = Randomizer.nextBoolean();
        
        // Get current (old) attachment height
        double oldHeight;
        if (moveDeparture) {
            oldHeight = conv.getHeight1();
        } else {
            oldHeight = conv.getHeight2();
        }
        
        // Choose window:
        double w = apertureSizeInput.get()*acg.getRoot().getHeight();

        // Set new height
        double newHeight = oldHeight + (Randomizer.nextDouble() - 0.5)*w;
        
        // Check for boundary violation
        if (moveDeparture) {
            if (newHeight>conv.getHeight2() || newHeight>acg.getRoot().getHeight())
                return Double.NEGATIVE_INFINITY;
        } else {
            if (newHeight<conv.getHeight1())
                return Double.NEGATIVE_INFINITY;
        }
        
        // Get node below current (old) attachment point
        Node nodeBelow;
        if (moveDeparture)
            nodeBelow = conv.getNode1();
        else
            nodeBelow = conv.getNode2();

        // Choose node below new attachment point
        if (newHeight<oldHeight) {
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
            conv.setHeight1(newHeight);
            conv.setNode1(nodeBelow);
        } else {
            conv.setHeight2(newHeight);
            conv.setNode2(nodeBelow);
        }
        
        return logHR;
    }
    
}
