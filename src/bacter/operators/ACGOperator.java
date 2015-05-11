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
import beast.core.Input;
import beast.core.Operator;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import feast.input.In;

/**
 * Abstract class of operators which act on the ConversionGraph state.
 * 
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public abstract class ACGOperator extends Operator {
    
    public Input<ConversionGraph> argInput = new In<ConversionGraph>(
            "acg", "Ancestral conversion graph.").setRequired();
    
    protected ConversionGraph acg;

    protected class ProposalFailed extends Exception {
        public ProposalFailed() { }
    }

    @Override
    public void initAndValidate() throws Exception {
        acg = argInput.get();
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
     * @return true if region map is valid under the restricted CO model.
     */
    protected boolean regionMapIsValid() {
        for (Alignment alignment : acg.getAlignments()) {
            for (int c = 0; c < acg.getConversions(alignment).size() - 1; c++) {
                Conversion thisConv = acg.getConversions(alignment).get(c);
                Conversion nextConv = acg.getConversions(alignment).get(c + 1);
                if (nextConv.getStartSite() - thisConv.getEndSite() < 2)
                    return false;
            }
        }
        return true;
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

        for (Alignment alignment : acg.getAlignments()) {
            for (Conversion conv : acg.getConversions(alignment)) {
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

        for (Alignment alignment : acg.getAlignments()) {
            for (Conversion conv : acg.getConversions(alignment)) {
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
        for (Alignment alignment : acg.getAlignments()) {
            if (idx<acg.getConvCount(alignment))
                return acg.getConversions(alignment).get(idx);
            else
                idx -= acg.getConvCount(alignment);
        }

        throw new IllegalStateException("Programmer error: loop fell through" +
                " in chooseConversion().");
    }

    /**
     * @return alignment selected randomly proportional to its length
     */
    protected Alignment chooseAlignment() {

        int z = Randomizer.nextInt(acg.getTotalSequenceLength());
        Alignment alignment = null;
        for (Alignment potentialAlignment : acg.getAlignments()) {
            if (z < potentialAlignment.getSiteCount())
                return alignment;
            else
                z -= potentialAlignment.getSiteCount();
        }

        throw new IllegalStateException("Programmer error: loop fell through" +
                " in chooseAlignment().");
    }
}
