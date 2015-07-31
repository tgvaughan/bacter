/*
 * Copyright (C) 2015 Tim Vaughan (tgvaughan@gmail.com)
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

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
@Description("Implementation of subtree slide that aims to sensibly deal " +
        "with conversions")
public class CFSubtreeSlide extends CFOperator {

    public Input<Double> sizeInput = new Input<>("size",
            "Size of window slide is confined to.", 1.0);

    public Input<Boolean> useGaussianInput = new Input<>("useGaussian",
            "Whether to use gaussian (true, default) or uniform delta.", true);

    @Override
    public double proposal() {

        double logHGF = 0.0;
        double logHalf = Math.log(0.5);

        // Choose non-root node:
        Node node = acg.getNode(Randomizer.nextInt(acg.getNodeCount()-1));

        Node parent = node.getParent();
        Node sister = getSibling(node);
        
        double delta = useGaussianInput.get()
                ? Randomizer.nextGaussian()*sizeInput.get()
                : (Randomizer.nextDouble()-0.5)*sizeInput.get();

        double newHeight = parent.getHeight() + delta;

        // Reject invalid moves:
        if (newHeight<node.getHeight())
            return Double.NEGATIVE_INFINITY;

        if (delta<0) {

            Node newSister = sister;
            while (newHeight<newSister.getHeight()) {
                if (newSister.isLeaf())
                    return Double.NEGATIVE_INFINITY;

                newSister = Randomizer.nextBoolean()
                        ? newSister.getLeft()
                        : newSister.getRight();

                logHGF -= logHalf;
            }

            // Update topology:
            logHGF += collapseConversions(node, newSister, newHeight);

        } else {

            // BUG HERE: need to skip nodeParent
            Node newSister = sister;
            while (!newSister.isRoot()
                        && newHeight>newSister.getParent().getHeight()) {
                    newSister = newSister.getParent();
                    logHGF += logHalf;
            }

            // Topology modification
            logHGF -= expandConversions(node, newSister, newHeight);

        }

        if (acg.isInvalid())
            throw new RuntimeException("CFCS proposed invalid state.");

        return logHGF;
    }
    
}
