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

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;

import java.util.Random;

/**
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
@Description("Implementation of subtree slide that aims to sensibly deal " +
        "with conversions")
public class CFSubtreeSlide extends CFOperator {

    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "Node height will be scaled by a factor between [scaleFactor,1/scaleFactor].",
            0.8);

    @Override
    public double proposal() {

        double logHGF = 0.0;
        double logHalf = Math.log(0.5);

        // Choose non-root node:
        Node srcNode = acg.getNode(Randomizer.nextInt(acg.getNodeCount() - 1));

        Node srcNodeP = srcNode.getParent();
        Node srcNodeS = getSibling(srcNode);

        double fMin = Math.min(scaleFactorInput.get(), 1.0 / scaleFactorInput.get());
        double f = fMin + Randomizer.nextDouble()*(1.0/fMin - fMin);
        logHGF += -Math.log(f);

        double newHeight = srcNodeP.getHeight()*f;

        // Reject invalid moves:
        if (newHeight<srcNode.getHeight())
            return Double.NEGATIVE_INFINITY;

        if (f<1.0) {

            // Search downwards for new attachment point:
            Node newSister = srcNodeS;
            while (newHeight<newSister.getHeight()) {
                if (newSister.isLeaf())
                    return Double.NEGATIVE_INFINITY;

                newSister = Randomizer.nextBoolean()
                        ? newSister.getLeft()
                        : newSister.getRight();

                logHGF -= logHalf;
            }

            // Update topology:
            logHGF += collapseConversions(srcNode, newSister, newHeight);

        } else {

            // Search upwards for new attachment point:
            Node newSister = srcNodeS;
            while (getParentSkip(newSister, srcNodeP) != null
                        && newHeight> getParentSkip(newSister, srcNodeP).getHeight()) {
                    newSister = getParentSkip(newSister, srcNodeP);
                    logHGF += logHalf;
            }

            // Update topology
            logHGF -= expandConversions(srcNode, newSister, newHeight);

        }

        assert !acg.isInvalid() : "CFSTS proposed invalid state.";

        return logHGF;
    }

    /**
     * Return parent of node unless parent==skip in which case grandparent
     * is returned.
     *
     * @param node node to find parent of
     * @param skip parent node to skip
     * @return conditioned parent
     */
    private Node getParentSkip(Node node, Node skip) {
        Node parent = node.getParent();
        if (parent == skip)
            return parent.getParent();
        else
            return parent;
    }
    
}
