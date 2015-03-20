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
import bacter.operators.EdgeCreationOperator;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.math.GammaFunction;
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

    //int count = 0;

    @Override
    public double proposal() {

        //System.out.println(++count);
        //System.out.println(acg.getExtendedNewick(true));

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

            System.out.println(acg.getExtendedNewick(true));

            double logHGF = 0.0;

            double t_srcNodeG = srcNodeP.getParent().getHeight();

            logHGF += Math.log(1.0/(t_srcNodeG - Math.max(t_srcNode, t_srcNodeS)));

            double newTime = t_destNode
                    + Randomizer.nextExponential(1.0/(alpha*t_destNode));

            logHGF -= Math.log(1.0/(alpha*t_destNode))
                    - (1.0/alpha)*(newTime/t_destNode - 1.0);

            disconnectEdge(srcNode);
            connectEdge(srcNode, destNode, newTime);
            acg.setRoot(srcNodeP);

            // TODO: Randomly reconnect some of the conversions ancestral
            // to srcNode to the new edge above srcNode.

            // Create conversions between the edge above destNode and the
            // srcNode above t_destNode.
            double L = 2*(newTime - t_destNode);
            double Nexp = L*rhoInput.get().getValue()
                    *(acg.getSequenceLength() + deltaInput.get().getValue());
            int N = (int)Randomizer.nextPoisson(Nexp);

            logHGF -= -Nexp + N*Math.log(Nexp);
                    //- GammaFunction.lnGamma(N+1);

            for (int i=0; i<N; i++) {
                Conversion conv = new Conversion();
                double u = L*Randomizer.nextDouble();
                if (u < 0.5*L) {
                    conv.setNode1(destNode);
                    conv.setHeight1(t_destNode + u);
                } else {
                    conv.setNode1(srcNode);
                    conv.setHeight2(t_destNode + (u - 0.5*L));
                }
                logHGF -= Math.log(1.0/L) + coalesceEdge(conv) + drawAffectedRegion(conv);
                acg.addConversion(conv);
            }

            // DEBUG
            if (!acg.isValid())
                return Double.NEGATIVE_INFINITY;

            System.out.println(acg.getExtendedNewick(true));
            return logHGF;
        }

        if (srcNodeP.isRoot()) {
            // Backward root move

            if (!includeRootInput.get())
                return Double.NEGATIVE_INFINITY;

            double logHGF = 0.0;

            // Remove conversions that will become root loops on transformation
            List<Conversion> toRemove = new ArrayList<>();
            for (Conversion conv : acg.getConversions()) {
                if ((conv.getNode1() == srcNode || conv.getNode1() == srcNodeS)
                        && conv.getHeight1() > t_srcNodeS) {
                    toRemove.add(conv);
                }
            }

            double L = 2*(t_srcNodeP - t_srcNodeS);

            double Nexp = L*rhoInput.get().getValue()
                    *(acg.getSequenceLength()+deltaInput.get().getValue());
            logHGF += -Nexp + toRemove.size()*Math.log(Nexp);
                    //- GammaFunction.lnGamma(toRemove.size()+1);

            for (Conversion conv : toRemove) {
                logHGF += Math.log(1.0/L) + getEdgeCoalescenceProb(conv) + getAffectedRegionProb(conv);
                acg.deleteConversion(conv);
            }


            logHGF += Math.log(1.0/(alpha*t_srcNodeS))
                    - (1.0/alpha)*(t_srcNodeP/t_srcNodeS - 1.0);

            double min_newTime = Math.max(t_srcNode, t_destNode);
            double t_destNodeP = destNodeP.getHeight();
            double newTime = min_newTime
                    + (t_destNodeP - min_newTime)*Randomizer.nextDouble();

            logHGF -= Math.log(1.0/(t_destNodeP - min_newTime));

            disconnectEdge(srcNode);
            connectEdge(srcNode, destNode, newTime);
            acg.setRoot(srcNodeS);

            // TODO: Reconnect conversions on edge above srcNode older than
            // newTime to edges ancestral to destNode.

            // DEBUG
            if (!acg.isValid())
                return Double.NEGATIVE_INFINITY;

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

        disconnectEdge(srcNode);
        connectEdge(srcNode, destNode, newTime);

        // DEBUG
        if (!acg.isValid())
            return Double.NEGATIVE_INFINITY;

        return logHGF;
    }
    
    /**
     * Returns true if srcNode CANNOT be used for the WB move.
     *
     * @param srcNode
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
     * @param srcNode
     * @param destNode
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
