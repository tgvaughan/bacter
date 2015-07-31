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
import bacter.ConversionGraph;
import bacter.Locus;
import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * Abstract class of operators which act on the ConversionGraph state.
 * 
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public abstract class ACGOperator extends Operator {
    
    public Input<ConversionGraph> acgInput = new Input<>(
            "acg",
            "Ancestral conversion graph.",
            Input.Validate.REQUIRED);
    
    protected ConversionGraph acg;

    @Override
    public void initAndValidate() throws Exception {
        acg = acgInput.get();
    }
    
    /**
     * Return sister of node.
     * 
     * @param node to return sister of
     * @return sister of node
     */
    protected Node getSibling(Node node) {
        Node parent = node.getParent();
        if (parent.getLeft() == node)
            return parent.getRight();
        else
            return parent.getLeft();
    }

    /**
     * Disconnect edge <node, node.parent> from its sister and
     * grandparent (if the grandparent exists), forming a new edge
     * <sister, grandparent>. All conversions originally above node.parent
     * are re-attached to sister.
     *
     * Conversions on edge above node are not modified.
     * 
     * @param node base of edge to detach.
     */
    protected void disconnectEdge(Node node) {

        if (node.isRoot())
            throw new IllegalArgumentException("Programmer error: "
                    + "root argument passed to disconnectEdge().");

        Node parent = node.getParent();
        Node sister = getSibling(node);

        if (parent.isRoot()) {
            parent.removeChild(sister);
            sister.setParent(null);
        } else {
            Node grandParent = parent.getParent();
            grandParent.removeChild(parent);
            parent.setParent(null);
            parent.removeChild(sister);
            grandParent.addChild(sister);
        }

        for (Locus locus : acg.getLoci()) {
            for (Conversion conv : acg.getConversions(locus)) {
                if (conv.getNode1() == parent)
                    conv.setNode1(sister);

                if (conv.getNode2() == parent)
                    conv.setNode2(sister);
            }
        }
    }

    /**
     * Connect edge <node, node.parent> above destEdgeBase, forming new
     * edge <destEdgeBase, node.parent> and <node.parent, destEdgeBase.parent>.
     * All conversions above destEdgeBase that are older than destTime
     * are transferred to the new edge above node.parent.
     *
     * Conversions on edge above node are not modified.
     * 
     * @param node base of edge to attach
     * @param destEdgeBase base of edge to be bisected
     * @param destTime height at which bisection will occur
     */
    protected void connectEdge(Node node, Node destEdgeBase, double destTime) {

        if (node.isRoot())
            throw new IllegalArgumentException("Programmer error: "
                    + "root argument passed to connectEdge().");

        Node parent = node.getParent();
        
        if (destEdgeBase.isRoot()) {
            parent.addChild(destEdgeBase);
        } else {
            Node grandParent = destEdgeBase.getParent();
            grandParent.removeChild(destEdgeBase);
            grandParent.addChild(parent);
            parent.addChild(destEdgeBase);
        }

        parent.setHeight(destTime);

        for (Locus locus : acg.getLoci()) {
            for (Conversion conv : acg.getConversions(locus)) {
                if (conv.getNode1() == destEdgeBase && conv.getHeight1() > destTime)
                    conv.setNode1(parent);

                if (conv.getNode2() == destEdgeBase && conv.getHeight2() > destTime)
                    conv.setNode2(parent);
            }
        }
    }

    /**
     * @return conversion selected uniformly at random
     */
    protected Conversion chooseConversion() {
        int idx = Randomizer.nextInt(acg.getTotalConvCount());
        for (Locus locus : acg.getLoci()) {
            if (idx<acg.getConvCount(locus))
                return acg.getConversions(locus).get(idx);
            else
                idx -= acg.getConvCount(locus);
        }

        throw new IllegalStateException("Programmer error: loop fell through" +
                " in chooseConversion().");
    }

    /**
     * @return alignment selected proportional to its length
     */
    protected Locus chooseLocus() {
        int z = Randomizer.nextInt(acg.getTotalSequenceLength());
        for (Locus locus : acg.getLoci()) {
            if (z < locus.getSiteCount())
                return locus;
            else
                z -= locus.getSiteCount();
        }

        throw new IllegalStateException("Programmer error: loop fell through" +
                " in chooseAlignment().");
    }

    /**
     * @return maximum height of children of root
     */
    protected double getMaxRootChildHeight() {
        double max = 0.0;
        for (Node child : acg.getRoot().getChildren())
            max = Math.max(max, child.getHeight());

        return max;
    }
}
