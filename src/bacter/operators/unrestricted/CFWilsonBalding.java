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
package bacter.operators.unrestricted;

import bacter.Conversion;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import feast.input.In;

import java.util.ArrayList;
import java.util.List;

/**
 * Implementation of Wilson-Balding operator modified for the clonal frame
 * of the ACG.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Wilson-Balding operator for ACG clonal frames.")
public class CFWilsonBalding extends ConversionCreationOperator {

    public Input<Double> alphaInput = new In<Double>("alpha", "Root height "
            + "proposal parameter").setRequired();

    public Input<Boolean> includeRootInput = new In<Boolean>("includeRoot",
            "Whether to include root variants of move.").setDefault(true);

    public Input<RealParameter> rhoInput = new In<RealParameter>("rho",
            "Conversion rate.").setRequired();

    private double alpha;

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        alpha = alphaInput.get();
    }

//    int count = 0;

    @Override
    public double proposal() {

//        System.out.println(++count);

//        System.out.println(acg.getExtendedNewick(true));

        // Determine whether we can apply this operator:
        if (acg.getNodeCount()<3)
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

            // DEBUG
            if (!acg.isValid())
                throw new IllegalStateException("CFWB proposed invalid state.");

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

            // DEBUG
            if (!acg.isValid())
                throw new IllegalStateException("CFWB proposed invalid state.");

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

        // DEBUG
        if (!acg.isValid())
            throw new IllegalStateException("CFWB proposed invalid state.");

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

    /**
     * Take conversions which connect to edge above srcNode at times greater than
     * destTime and attach them instead to the lineage above destNode.
     *
     * Assumes topology has not yet been altered.
     *
     * @param srcNode   source node for move
     * @param destNode  dest node for move
     * @param destTime  new time of attachment of edge above srcNode to edge
     *                  above destNode
     * @return log probability of the collapsed attachments.
     */
    private double collapseConversions(Node srcNode, Node destNode, double destTime) {
        double logP = 0.0;

        if (destNode.isRoot())
            throw new IllegalArgumentException("Destination node should " +
                    "never be root in argument to collapseConversions.");

        boolean reverseRootMove = srcNode.getParent().isRoot();
        Node srcNodeS = getSibling(srcNode);

        // Collapse non-root conversions

        Node node = destNode;
        while (!node.isRoot() &&
                node.getHeight() < srcNode.getParent().getHeight()) {
            for (Conversion conv : acg.getConversions()) {
                if (conv.getNode1() == srcNode
                        && conv.getHeight1() > Math.max(destTime, node.getHeight())
                        && conv.getHeight1() < node.getParent().getHeight()) {
                    conv.setNode1(node);
                    logP += Math.log(0.5);
                }

                if (conv.getNode2() == srcNode
                        && conv.getHeight2() > Math.max(destTime, node.getHeight())
                        && conv.getHeight2() < node.getParent().getHeight()) {
                    conv.setNode2(node);
                    logP += Math.log(0.5);
                }
            }

            node = node.getParent();
        }

        // Collapse conversions between srcNode edge and its sibling
        // if this was a reverse root move.

        if (reverseRootMove) {
            double maxChildHeight = Math.max(srcNodeS.getHeight(), srcNode.getHeight());

            double L = 2.0*(srcNode.getParent().getHeight() - maxChildHeight);

            double Nexp = L*rhoInput.get().getValue()
                    *(acg.getSequenceLength() + deltaInput.get().getValue());

            List<Conversion> toRemove = new ArrayList<>();
            for (Conversion conv : acg.getConversions()) {
                if (conv.getNode1() == srcNodeS
                        || (conv.getNode1() == srcNode && conv.getHeight1()>maxChildHeight))
                    toRemove.add(conv);
            }

            logP += -Nexp + toRemove.size()*Math.log(Nexp);
            // Factorial cancelled due to sum over permutations of
            // individual conversion states

            for (Conversion conv : toRemove) {
                logP += Math.log(1.0/L)
                        + getAffectedRegionProb(conv)
                        + getEdgeCoalescenceProb(conv);

                acg.deleteConversion(conv);
            }
        }

        // Apply topology modifications.
        disconnectEdge(srcNode);
        connectEdge(srcNode, destNode, destTime);

        if (reverseRootMove)
            acg.setRoot(srcNodeS);

        return logP;
    }

    /**
     * Take length of new edge above srcNode that is greater than the
     * original height of srcNode.parent and shifts a random fraction of
     * conversion attachments to it from the lineage above destNode.
     *
     * In the case that destNode was the root, the conversions starting
     * above destNode are drawn from the prior.
     *
     * Assumes topology has not yet been altered.
     *
     * @param srcNode source node for the move
     * @param destNode dest node for the move
     * @param destTime new time drawn for srcNode.P.
     * @return log probability of new conversion configuration.
     */
    private double expandConversions(Node srcNode, Node destNode, double destTime) {
        double logP = 0.0;

        boolean forwardRootMove = destNode.isRoot();

        Node node = srcNode.getParent();
        while (!node.isRoot()) {
            for (Conversion conv : acg.getConversions()) {
                if (conv.getNode1() == node && conv.getHeight1() < destTime) {
                    if (Randomizer.nextBoolean())
                        conv.setNode1(srcNode);
                    logP += Math.log(0.5);
                }

                if (conv.getNode2() == node && conv.getHeight2() < destTime) {
                    if (Randomizer.nextBoolean())
                        conv.setNode2(srcNode);
                    logP += Math.log(0.5);
                }

            }

            node = node.getParent();
        }

        // Apply topology modifications.
        disconnectEdge(srcNode);
        connectEdge(srcNode, destNode, destTime);

        // Add conversions between srcNode edge and its sibling if
        // this was a forward root move
        if (forwardRootMove) {
            acg.setRoot(srcNode.getParent());

            double L = 2.0*(destTime - destNode.getHeight());
            double Nexp = L*rhoInput.get().getValue()
                    *(acg.getSequenceLength() + deltaInput.get().getValue());
            int N = (int)Randomizer.nextPoisson(Nexp);

            logP += -Nexp + N*Math.log(Nexp); // Factorial cancels

            for (int i=0; i<N; i++) {
                Conversion conv = new Conversion();

                double u = Randomizer.nextDouble()*L;
                if (u < 0.5*L) {
                    conv.setNode1(destNode);
                    conv.setHeight1(destNode.getHeight() + u);
                } else {
                    conv.setNode1(srcNode);
                    conv.setHeight1(destNode.getHeight() + u - 0.5*L);
                }

                logP += Math.log(1.0/L) + drawAffectedRegion(conv) + coalesceEdge(conv);

                acg.addConversion(conv);
            }

        }

        return logP;
    }

}
