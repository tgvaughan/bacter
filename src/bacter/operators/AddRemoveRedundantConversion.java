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
@Description("Operator which adds and removes redundant conversions to/from an ACG, " +
        "i.e. those which do not alter the CF topology.")
public class AddRemoveRedundantConversion extends ACGOperator {

    public Input<Double> apertureInput = new Input<>("aperture",
            "Maximum aperture used to select conversion attachment points, " +
                    "expressed as a fraction of node age.", 0.1);

    public AddRemoveRedundantConversion() { }
    
    @Override
    public double proposal() {
        double logHGF = 0;

        // Select non-root CF node
        Node cfNode = acg.getNode(Randomizer.nextInt(acg.getNodeCount()-1));
        Node cfParent = cfNode.getParent();

        double maxL = apertureInput.get()*acg.getRoot().getHeight();

        if (Randomizer.nextBoolean()) {
            
            // Add redundant conversion

            Conversion newConv = new Conversion();

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
                    if (fromPoint < Math.min(getEdgeLength(child), maxL)) {
                        newConv.setNode1(child);
                        newConv.setHeight1(cfNode.getHeight() - fromPoint);
                        break;
                    }
                    fromPoint -= Math.min(getEdgeLength(child), maxL);
                }

            }

            L = Math.min(getEdgeLength(cfParent), maxL);
            for (Node child : cfParent.getChildren())
                L += Math.min(getEdgeLength(child), maxL);

            logHGF -= Math.log(1.0/L);

            double toPoint = L*Randomizer.nextDouble();
            if (toPoint < Math.min(getEdgeLength(cfParent), maxL)) {
                newConv.setNode2(cfParent);
                newConv.setHeight2(cfParent.getHeight() + toPoint);
            } else {
                toPoint -= Math.min(getEdgeLength(cfParent), maxL);
                for (Node child : cfParent.getChildren()) {
                    if (toPoint < Math.min(getEdgeLength(child), maxL)) {
                        newConv.setNode2(child);
                        newConv.setHeight2(cfParent.getHeight() - toPoint);
                        break;
                    }
                    toPoint -= Math.min(getEdgeLength(child), maxL);
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

            double L = Math.min(getEdgeLength(cfNode), maxL);
            for (Node child : cfNode.getChildren())
                L += Math.min(getEdgeLength(child), maxL);
            logHGF += Math.log(1.0/L);

            L = Math.min(getEdgeLength(cfParent), maxL);
            for (Node child : cfParent.getChildren())
                    L += Math.min(getEdgeLength(child), maxL);
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
        Node cfParent = cfNode.getParent();

        List<Conversion> redundantConversions = new ArrayList<>();

        double maxL = acg.getRoot().getHeight()*apertureInput.get();

        for (Locus locus : acg.getLoci()) {
            for (Conversion conv : acg.getConversions(locus)) {

                if (((conv.getNode1() == cfNode || conv.getNode1().getParent() == cfNode)
                        && Math.abs(conv.getHeight1()-cfNode.getHeight()) < maxL)
                    && (conv.getNode2() == cfParent || conv.getNode2().getParent() == cfParent)
                        && Math.abs(conv.getHeight2()-cfParent.getHeight()) < maxL) {

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

    private double drawAffectedRegion(Conversion conv) {
        double logP = 0.0;

        Locus locus = acg.getLoci().get(Randomizer.nextInt(acg.getLoci().size()));
        conv.setLocus(locus);
        logP += Math.log(1.0/acg.getLoci().size());

        if (!acg.wholeLocusModeOn()) {
            int site1 = Randomizer.nextInt(locus.getSiteCount());
            int site2 = Randomizer.nextInt(locus.getSiteCount());

            if (site1 < site2) {
                conv.setStartSite(site1);
                conv.setEndSite(site2);
            } else {
                conv.setStartSite(site2);
                conv.setEndSite(site1);
            }

            logP += 2.0 * Math.log(1.0 / locus.getSiteCount());
            if (site1 != site2)
                logP += Math.log(2.0);
        } else {
            conv.setStartSite(0);
            conv.setEndSite(locus.getSiteCount()-1);
        }

        return logP;
    }

    private double getAffectedRegionProb(Conversion conv) {
        double logP = 0.0;

        logP += Math.log(1.0/acg.getLoci().size());

        if (!acg.wholeLocusModeOn()) {
            logP += 2.0 * Math.log(1.0 / conv.getLocus().getSiteCount());
            if (conv.getStartSite() != conv.getEndSite())
                logP += Math.log(2.0);
        }

        return logP;
    }
}
