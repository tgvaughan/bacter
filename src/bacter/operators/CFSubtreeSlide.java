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
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * Implementation of subtree slide that aims to sensibly deal with
 * conversions.
 *
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
public class CFSubtreeSlide extends ACGOperator {

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
            while (true) {
                for (Conversion conv : acg.getConversions()) {
                    if (conv.getNode1() == node
                            && conv.getHeight1()>newHeight
                            && conv.getHeight1()>newSister.getHeight()) {
                        conv.setNode1(newSister);
                        logHGF += logHalf;
                    }
                    
                    if (conv.getNode2() == node
                            && conv.getHeight2()>newHeight
                            && conv.getHeight2()>newSister.getHeight()) {
                        conv.setNode2(newSister);
                        logHGF += logHalf;
                    }
                }
                
                if (newHeight<newSister.getHeight()) {
                    if (newSister.isLeaf())
                        return Double.NEGATIVE_INFINITY;
                    
                    newSister = Randomizer.nextBoolean()
                            ? newSister.getLeft()
                            : newSister.getRight();

                    logHGF -= logHalf;
                } else
                    break;
            }

            // Update topology:
            if (newSister != sister) {
                if (parent.isRoot()) {
                    parent.removeChild(sister);
                    sister.setParent(null);
                    acg.setRoot(sister);
                } else {
                    Node grandParent = parent.getParent();
                    grandParent.removeChild(parent);
                    parent.removeChild(sister);
                    grandParent.addChild(sister);
                }

                Node grandParent = newSister.getParent();
                grandParent.removeChild(newSister);
                grandParent.addChild(parent);
                parent.addChild(newSister);
            }

            parent.setHeight(newHeight);

            for (Conversion conv : acg.getConversions()) {
                if (conv.getNode1() == newSister && conv.getHeight1()>newHeight) {
                    if (!parent.isRoot())
                        conv.setNode1(parent);
                    else
                        return Double.NEGATIVE_INFINITY;
                }

                if (conv.getNode2() == newSister && conv.getHeight2()>newHeight)
                    conv.setNode2(parent);
            }
                
        } else {

            System.out.println(acg.getExtendedNewick(true));
            
            // Disconnect <node,parent>
            if (parent.isRoot()) {
                parent.removeChild(sister);
                sister.setParent(null);
                acg.setRoot(sister);
            } else {
                Node grandParent = parent.getParent();
                grandParent.removeChild(parent);
                parent.setParent(null);
                parent.removeChild(sister);
                grandParent.addChild(sister);
            }

            // Move any conversions on edge previously above parent to edge
            // above sister
            for (Conversion conv: acg.getConversions()) {
                if (conv.getNode1() == parent)
                    conv.setNode1(sister);
                if (conv.getNode2() == parent)
                    conv.setNode2(sister);
            }

            Node newSister = sister;
            while (true) {
                for (Conversion conv : acg.getConversions()) {
                    if (conv.getNode1() == newSister
                            && conv.getHeight1()>parent.getHeight()
                            && conv.getHeight1()<newHeight) {
                        if (Randomizer.nextBoolean())
                            conv.setNode1(node);
                        
                        logHGF -= logHalf;
                    }
                    
                    if (conv.getNode2() == newSister
                            && conv.getHeight2()>parent.getHeight()
                            && conv.getHeight2()<newHeight) {
                        if (Randomizer.nextBoolean())
                            conv.setNode2(node);
                        
                        logHGF -= logHalf;
                    }
                }
                
                if (!newSister.isRoot()
                        && newHeight>newSister.getParent().getHeight()) {
                    newSister = newSister.getParent();
                    logHGF += logHalf;
                } else
                    break;
            }

            // Topology modification
            if (newSister.isRoot()) {
                parent.addChild(newSister);
                acg.setRoot(parent);
            } else {
                Node grandParent = newSister.getParent();
                grandParent.removeChild(newSister);
                grandParent.addChild(parent);
                parent.addChild(newSister);
            }

            parent.setHeight(newHeight);

            for (Conversion conv : acg.getConversions()) {
                if (conv.getNode1()==newSister && conv.getHeight2()<newHeight)
                    conv.setNode1(parent);

                if (conv.getNode2()==newSister && conv.getHeight2()<newHeight)
                    conv.setNode2(parent);
            }

        }

        if (delta<0)
            return Double.NEGATIVE_INFINITY;

        System.out.println(acg.getExtendedNewick(true));
        return logHGF;
    }
    
}
