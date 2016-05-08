/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
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
 * Implementation of Wilson-Balding operator modified for the clonal frame
 * of the ACG.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Wilson-Balding operator for ACG clonal frames.")
public class CFWilsonBalding extends CFOperator {

    public Input<Double> alphaInput = new Input<>("alpha", "Root height "
            + "proposal parameter", Input.Validate.REQUIRED);

    public Input<Boolean> includeRootInput = new Input<>("includeRoot",
            "Whether to include root variants of move.", true);

    private double alpha;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        alpha = alphaInput.get();
    }

    @Override
    public double proposal() {

        // Determine whether we can apply this operator:
        if (acg.getLeafNodeCount()<3)
            return Double.NEGATIVE_INFINITY;

        // Select non-root node:
        Node srcNode;
        do {
            srcNode = acg.getNode(Randomizer.nextInt(acg.getNodeCount()-1));
        } while (invalidSrcNode(srcNode));
        Node srcNodeP = srcNode.getParent();
        Node srcNodeS = getSibling(srcNode);
        double t_srcNode = srcNode.getHeight();
        double t_srcNodeP = srcNodeP.getHeight();
        double t_srcNodeS = srcNodeS.getHeight();

        // Select destination branch node:
        Node destNode;
        do {
            destNode = acg.getNode(Randomizer.nextInt(acg.getNodeCount()));
        } while (invalidDestNode(srcNode, destNode));
        Node destNodeP = destNode.getParent();
        double t_destNode = destNode.getHeight();

        if (destNode.isRoot()) {
            // Forward root move

            if (!includeRootInput.get())
                return Double.NEGATIVE_INFINITY;

            double logHGF = 0.0;

            double t_srcNodeG = srcNodeP.getParent().getHeight();

            logHGF += Math.log(1.0/(t_srcNodeG - Math.max(t_srcNode, t_srcNodeS)));

            double newTime = t_destNode
                    + Randomizer.nextExponential(1.0/(alpha*t_destNode));

            logHGF -= Math.log(1.0/(alpha*t_destNode))
                    - (1.0/alpha)*(newTime/t_destNode - 1.0);

            // Randomly reconnect some of the conversions ancestral
            // to srcNode to the new edge above srcNode.
            logHGF -= expandConversions(srcNode, destNode, newTime);

            assert !acg.isInvalid() : "CFWB proposed invalid state.";

            return logHGF;
        }

        if (srcNodeP.isRoot()) {
            // Backward root move

            if (!includeRootInput.get())
                return Double.NEGATIVE_INFINITY;

            double logHGF = 0.0;

            logHGF += Math.log(1.0/(alpha*t_srcNodeS))
                    - (1.0/alpha)*(t_srcNodeP/t_srcNodeS - 1.0);

            double min_newTime = Math.max(t_srcNode, t_destNode);
            double t_destNodeP = destNodeP.getHeight();
            double newTime = min_newTime
                    + (t_destNodeP - min_newTime)*Randomizer.nextDouble();

            logHGF -= Math.log(1.0/(t_destNodeP - min_newTime));

            // Reconnect conversions on edge above srcNode older than
            // newTime to edges ancestral to destNode.
            logHGF += collapseConversions(srcNode, destNode, newTime);

            assert !acg.isInvalid() : "CFWB proposed invalid state.";

            return logHGF;
        }

        // Non-root move

        double logHGF = 0.0;

        double t_srcNodeG = srcNodeP.getParent().getHeight();

        logHGF += Math.log(1.0/(t_srcNodeG - Math.max(t_srcNode, t_srcNodeS)));

        double min_newTime = Math.max(t_destNode, t_srcNode);
        double t_destNodeP = destNodeP.getHeight();
        double newTime = min_newTime
                + (t_destNodeP - min_newTime)*Randomizer.nextDouble();

        logHGF -= Math.log(1.0/(t_destNodeP - min_newTime));

        if (newTime < srcNodeP.getHeight())
            logHGF += collapseConversions(srcNode, destNode, newTime);
        else
            logHGF -= expandConversions(srcNode, destNode, newTime);

        assert !acg.isInvalid() : "CFWB proposed invalid state.";

        return logHGF;
    }
    
    /**
     * Returns true if srcNode CANNOT be used for the WB move.
     *
     * @param srcNode source node for move
     * @return True if srcNode invalid.
     */
    private boolean invalidSrcNode(Node srcNode) {

        if (srcNode.isRoot())
            return true;

        Node parent = srcNode.getParent();

        // This check is important for avoiding situations where it is
        // impossible to choose a valid destNode:
        if (parent.isRoot()) {

            Node sister = getSibling(srcNode);

            if (sister.isLeaf())
                return true;

            if (srcNode.getHeight() >= sister.getHeight())
                return true;
        }

        return false;
    }

    /**
     * Returns true if destNode CANNOT be used for the WB move in conjunction
     * with srcNode.
     *
     * @param srcNode   source node for move
     * @param destNode  destination node for move
     * @return True if destNode invalid.
     */
    private boolean invalidDestNode(Node srcNode, Node destNode) {

        if (destNode==srcNode
                || destNode==srcNode.getParent()
                || destNode.getParent()==srcNode.getParent())
            return true;

        Node destNodeP = destNode.getParent();

        return destNodeP != null && (destNodeP.getHeight() <= srcNode.getHeight());
    }
}
