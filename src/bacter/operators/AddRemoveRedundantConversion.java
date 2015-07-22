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

import bacter.Conversion;
import bacter.Locus;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which adds and removes redundant conversions to/from an ACG.")
public class AddRemoveRedundantConversion extends ConversionCreationOperator {

    public Input<Double> apertureInput = new Input<Double>("aperture",
            "Maximum aperture used to select conversion attachment points, " +
                    "expressed as a fraction of node age.", 0.1);

    public AddRemoveRedundantConversion() { }
    
    @Override
    public double proposal() {
        double logHGF = 0;

        // Select non-root CF node
        Node cfNode = acg.getNode(Randomizer.nextInt(acg.getNodeCount()-1));
        Node cfParent = cfNode.getParent();

        if (Randomizer.nextBoolean()) {
            
            // Add redundant conversion

            Conversion newConv = new Conversion();

            double maxL = apertureInput.get()*cfNode.getHeight();
            double L = Math.min(getEdgeLength(cfNode), maxL);
            for (Node child : cfNode.getChildren())
                L += Math.min(getEdgeLength(child), maxL);

            logHGF -= Math.log(1.0/L);

            double fromPoint = L*Randomizer.nextDouble();
            if (fromPoint < Math.min(getEdgeLength(cfNode), maxL)) {
                newConv.setNode1(cfNode);
                newConv.setHeight1(cfNode.getHeight() + fromPoint);
            } else {
                fromPoint -= Math.min(getEdgeLength(cfNode), maxL);
                for (Node child : cfNode.getChildren()) {
                    if (fromPoint < child.getLength()) {
                        newConv.setNode1(child);
                        newConv.setHeight1(child.getHeight() + fromPoint);
                        break;
                    }
                    fromPoint -= child.getLength();
                }

            }

            maxL = apertureInput.get()*cfParent.getHeight();
            for (Node child : cfParent.getChildren())
                L += Math.min(child.getLength(), maxL);

            logHGF -= Math.log(1.0/L);

            double toPoint = L*Randomizer.nextDouble();
            if (toPoint < cfParent.getLength()) {
                newConv.setNode2(cfParent);
                newConv.setHeight2(cfParent.getHeight() + toPoint);
            } else {
                toPoint -= cfParent.getLength();
                for (Node child : cfParent.getChildren()) {
                    if (toPoint < child.getLength()) {
                        newConv.setNode2(child);
                        newConv.setHeight2(child.getHeight() + toPoint);
                        break;
                    }
                    toPoint -= child.getLength();
                }

            }

            if (newConv.getHeight1()>newConv.getHeight2())
                return Double.NEGATIVE_INFINITY;

            logHGF -= drawAffectedRegion(newConv);

            // Add conversion
            acg.addConversion(newConv);

            // Add probability of reverse move deleting this conversion
            // to HGF:
            logHGF += Math.log(1.0/getRedundantConversions(cfNode).size());
        } else {
            
            // Remove
            
            // Identify redundant conversions
            List<Conversion> redundantConversions = getRedundantConversions(cfNode);

            if (redundantConversions.size() == 0)
                return Double.NEGATIVE_INFINITY;

            // Choose conversion to remove
            Conversion conv = redundantConversions.get(Randomizer.nextInt(redundantConversions.size()));
            logHGF -= Math.log(1.0/redundantConversions.size());

            // Add probability of reverse move generating this conversion
            // to HGF:

            double L = cfNode.getLength();
            for (Node child : cfNode.getChildren())
                L += child.getLength();
            logHGF += Math.log(1.0/L);

            L = cfParent.getLength();
            for (Node child : cfParent.getChildren())
                    L += child.getLength();
            logHGF += Math.log(1.0/L);

            logHGF += getAffectedRegionProb(conv);

            // Remove conversion
            acg.deleteConversion(conv);

        }

        return logHGF;
    }

    /**
     * Obtain list of redundant conversions.
     *
     * @param cfNode node identifying edge on CF
     * @return conversion list
     */
    private List<Conversion> getRedundantConversions(Node cfNode) {
        List<Conversion> redundantConversions = new ArrayList<>();

        for (Locus locus : acg.getLoci()) {
            for (Conversion conv : acg.getConversions(locus)) {
                if ((conv.getNode1() == cfNode || conv.getNode1().getParent() == cfNode) &&
                        (conv.getNode2() == cfNode.getParent() || (!conv.getNode2().isRoot()
                                && conv.getNode2().getParent() == cfNode.getParent()))) {
                    redundantConversions.add(conv);
                }
            }
        }

        return redundantConversions;
    }

    /**
     * Return length of edge above node.  Unlike node.getLength(), this method
     * sensibly handles the root.
     *
     * @param node node below edge
     * @return length of edge above node
     */
    private double getEdgeLength(Node node) {
        if (node.isRoot())
            return Double.POSITIVE_INFINITY;
        else
            return node.getParent().getHeight() - node.getHeight();
    }
}
