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

import bacter.Conversion;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * Uniform operator for clonal frame nodes. This operator is capable of
 * shifting an internal CF node past conversions, getting better acceptance
 * rates than the standard uniform operator when a large number of conversions
 * is present.
 *
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
public class CFUniform extends ACGOperator {

    @Override
    public double proposal() {

        double logHGF = 0.0;

        double logHalf = Math.log(0.5);

        // Select internal non-root node at random.
        Node node = acg.getNode(acg.getLeafNodeCount()
                + Randomizer.nextInt(acg.getInternalNodeCount()-1));

        Node parent = node.getParent();
        Node leftChild = node.getLeft();
        Node rightChild = node.getRight();

        // Choose new height uniformly between bounds:
        double minHeight = Math.max(leftChild.getHeight(), rightChild.getHeight());
        double maxHeight = parent.getHeight();
        double newHeight = minHeight + Randomizer.nextDouble()*(maxHeight-minHeight);

        if (newHeight>node.getHeight()) {
            for (Conversion conv : acg.getConversions()) {
                if (conv.getNode1() == node && conv.getHeight1()<newHeight) {
                    conv.setNode1(Randomizer.nextBoolean() ? leftChild : rightChild);
                    logHGF -= logHalf;
                }

                if (conv.getNode2() == node && conv.getHeight2()<newHeight) {
                    conv.setNode2(Randomizer.nextBoolean() ? leftChild : rightChild);
                    logHGF -= logHalf;
                }
            }
        } else {
            for (Conversion conv : acg.getConversions()) {
                if ((conv.getNode1() == leftChild || conv.getNode1() == rightChild)
                        && conv.getHeight1()>newHeight) {
                    conv.setNode1(node);
                    logHGF += logHalf;
                }

                if ((conv.getNode2() == leftChild || conv.getNode2() == rightChild)
                        && conv.getHeight2()>newHeight) {
                    conv.setNode2(node);
                    logHGF += logHalf;
                }
            }
        }

        node.setHeight(newHeight);

        return logHGF;
    }
    
}
