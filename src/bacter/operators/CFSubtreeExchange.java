package bacter.operators;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class CFSubtreeExchange extends CFOperator {

    public Input<Boolean> isNarrowInput = new Input<>("isNarrow",
            "Whether to use narrow exchange. (Default true.)", true);

    @Override
    public double proposal() {
        double logHGF = 0.0;

        Node srcNode, srcNodeP, destNode, destNodeP;
        if (isNarrowInput.get()) {

            // Narrow exchange selection:
            do {
                srcNode = acg.getNode(Randomizer.nextInt(acg.getNodeCount()-1));
            } while (srcNode.getParent().isRoot());
            srcNodeP = srcNode.getParent();
            destNode = getSibling(srcNodeP);
            destNodeP = destNode.getParent();

        } else {

            // Wide exchange selection:
            srcNode = acg.getNode(Randomizer.nextInt(acg.getNodeCount()-1));
            srcNodeP = srcNode.getParent();

            do {
                destNode = acg.getNode(Randomizer.nextInt(acg.getNodeCount()-1));
            } while(destNode == srcNode || destNode.getParent() == srcNode.getParent());
            destNodeP = destNode.getParent();
        }

        double t_srcNodeP = srcNodeP.getHeight();
        double t_destNodeP = destNodeP.getHeight();

        // Reject if substitution would result in negative branch lengths:
        if (destNode == srcNodeP || srcNode == destNodeP
                || destNode.getHeight()>t_srcNodeP
                || srcNode.getHeight()>t_destNodeP)
            return Double.NEGATIVE_INFINITY;

        if (t_srcNodeP > t_destNodeP) {
            logHGF += collapseConversions(srcNode, destNode, t_destNodeP);
            logHGF -= expandConversions(destNode, srcNode, t_srcNodeP);
        } else {
            logHGF -= expandConversions(srcNode, destNode, t_destNodeP);
            logHGF += collapseConversions(destNode, srcNode, t_srcNodeP);
        }

        assert !acg.isInvalid() : "CFSTX proposed invalid state.";

        return logHGF;
    }
}
