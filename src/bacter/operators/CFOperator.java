package bacter.operators;

import bacter.Conversion;
import bacter.Locus;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public abstract class CFOperator extends ConversionCreationOperator {

    public Input<RealParameter> rhoInput = new Input<>("rho",
            "Conversion rate.", Input.Validate.REQUIRED);

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
    protected double collapseConversions(Node srcNode, Node destNode, double destTime) {
        double logP = 0.0;

        boolean reverseRootMove = srcNode.getParent().isRoot();
        Node srcNodeP = srcNode.getParent();
        Node srcNodeS = getSibling(srcNode);
        double maxChildHeight = getMaxRootChildHeight();
        double volatileHeight = Math.max(maxChildHeight, destTime);

        // Collapse non-root conversions

        Node node = destNode;
        while (!node.isRoot() && node.getHeight() < srcNodeP.getHeight()) {

            double lowerBound = Math.max(destTime, node.getHeight());
            double upperBound = Math.min(node.getParent().getHeight(),
                    srcNodeP.getHeight());

            for (Locus locus : acg.getLoci()) {
                for (Conversion conv : acg.getConversions(locus)) {
                    if (conv.getHeight1() > lowerBound && conv.getHeight1() < upperBound) {
                        if (conv.getNode1() == srcNode)
                            conv.setNode1(node);

                        if (conv.getNode1() == node &&
                                (!reverseRootMove || conv.getHeight1() < volatileHeight))
                            logP += Math.log(0.5);
                    }

                    if (conv.getHeight2() > lowerBound && conv.getHeight2() < upperBound) {
                        if (conv.getNode2() == srcNode)
                            conv.setNode2(node);

                        if (conv.getNode2() == node
                                && (!reverseRootMove || conv.getNode1() != node
                                || conv.getHeight1() < volatileHeight))
                            logP += Math.log(0.5);
                    }
                }
            }

            node = node.getParent();
        }

        // Collapse conversions between srcNode edge and its sibling
        // if this was a reverse root move.

        if (reverseRootMove) {
            double L = 2.0*(srcNode.getParent().getHeight() - volatileHeight);

            double Nexp = L*rhoInput.get().getValue()
                    *(acg.getTotalSequenceLength()
                    + acg.getLoci().size()*(deltaInput.get().getValue() - 1.0));

            List<Conversion> toRemove = new ArrayList<>();
            for (Locus locus : acg.getLoci()) {
                for (Conversion conv : acg.getConversions(locus)) {
                    if (conv.getHeight1() > volatileHeight)
                        toRemove.add(conv);
                }
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

        if (reverseRootMove && destTime<maxChildHeight)
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
    protected double expandConversions(Node srcNode, Node destNode, double destTime) {
        double logP = 0.0;

        double volatileHeight = acg.getRoot().getHeight();
        boolean forwardRootMove = destTime > volatileHeight;

        Node node = srcNode.getParent();
        while (node != null) {
            for (Locus locus : acg.getLoci()) {
                for (Conversion conv : acg.getConversions(locus)) {
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

            double L = 2.0*(destTime - volatileHeight);
            double Nexp = L*rhoInput.get().getValue()
                    *(acg.getTotalSequenceLength()
                    + acg.getLoci().size()*(deltaInput.get().getValue() - 1.0));
            int N = (int)Randomizer.nextPoisson(Nexp);

            logP += -Nexp + N*Math.log(Nexp); // Factorial cancels

            for (int i=0; i<N; i++) {
                Conversion conv = new Conversion();

                double u = Randomizer.nextDouble()*L;
                if (u < 0.5*L) {
                    conv.setNode1(destNode);
                    conv.setHeight1(volatileHeight + u);
                } else {
                    conv.setNode1(srcNode);
                    conv.setHeight1(volatileHeight + u - 0.5*L);
                }

                logP += Math.log(1.0/L) + drawAffectedRegion(conv) + coalesceEdge(conv);

                acg.addConversion(conv);
            }

        }

        return logP;
    }

}
