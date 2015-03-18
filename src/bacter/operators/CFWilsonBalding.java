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

import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * Implementation of Wilson-Balding operator modified for the clonal frame
 * of the ACG.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class CFWilsonBalding extends ACGOperator {

    @Override
    public double proposal() {
        double logHGF = 0.0;

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

        }

        if (srcNodeP.isRoot()) {
            // Backward root move

        }

        // Non-root move

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
