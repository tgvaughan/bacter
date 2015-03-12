/*
 * Copyright (C) 2015 Tim Vaughan (tgvaughan@gmail.com)
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
import bacter.ConversionGraph;
import bacter.model.unrestricted.SimulatedACG;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.util.Randomizer;

/**
 * Operator which reversibly deletes a conversion, modifying the CF to match the
 * marginal tree of the original conversion.
 *
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
public class ClonalFrameConversionSwap extends ConversionCreationOperator {

    int count = 0;

    @Override
    public double proposal() {
        count += 1;

        if (Randomizer.nextBoolean()) {
            return deleteConversion();
        } else {
            return createConversion();
            //return Double.NEGATIVE_INFINITY;
        }
    }

    /**
     * Replaces a conversion with a modification to the clonal frame.
     * 
     * @return log Hastings-Green factor
     */
    public double deleteConversion() {
        double logHGF = 0.0;

        if (acg.getConvCount() == 0) {
            return Double.NEGATIVE_INFINITY;
        }

        Conversion chosenConv = acg.getConversions().get(
            Randomizer.nextInt(acg.getConvCount()));

        // Skip invisible conversions:
        if (chosenConv.getNode1() == chosenConv.getNode2()) {
            return Double.NEGATIVE_INFINITY;
        }

        logHGF -= Math.log(1.0 / acg.getConvCount());

        // Abort if conversions attach to node1 above height1.
        for (Conversion conv : acg.getConversions()) {
            if (conv == chosenConv) {
                continue;
            }

            if ((conv.getNode1() == chosenConv.getNode1()
                && conv.getHeight1() > chosenConv.getHeight1())
                || (conv.getNode2() == chosenConv.getNode1()
                && conv.getHeight2() > chosenConv.getHeight2())) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        // Move all conversion attachments from parent to sister.
        Node parent = chosenConv.getNode1().getParent();
        Node sister = getSibling(chosenConv.getNode1());

        for (Conversion conv : acg.getConversions()) {
            if (conv == chosenConv) {
                continue;
            }

            if (conv.getNode1() == parent) {
                conv.setNode1(sister);
            }

            if (conv.getNode2() == parent) {
                conv.setNode2(sister);
            }
        }

        if (chosenConv.getNode2() == parent) {
            chosenConv.setNode2(sister);
        }

        // Detach node1 from parent:
        parent.removeChild(sister);
        if (parent.getParent() != null) {
            Node grandParent = parent.getParent();
            grandParent.removeChild(parent);
            grandParent.addChild(sister);
            parent.setParent(null);
        } else {
            sister.setParent(null);
        }

        // Swap heights of parent and height2
        double oldParentHeight = parent.getHeight();
        parent.setHeight(chosenConv.getHeight2());
        chosenConv.setHeight2(oldParentHeight);

        // Move conversions from node2 to parent
        for (Conversion otherConv : acg.getConversions()) {
            if (otherConv == chosenConv) {
                continue;
            }

            if (otherConv.getNode1() == chosenConv.getNode2()
                && otherConv.getHeight1() > parent.getHeight()) {
                otherConv.setNode1(parent);
            }

            if (otherConv.getNode2() == chosenConv.getNode2()
                && otherConv.getHeight2() > parent.getHeight()) {
                otherConv.setNode2(parent);
            }
        }

        // Make final topology changes:
        Node oldParent = chosenConv.getNode2().getParent();
        parent.addChild(chosenConv.getNode2());
        if (oldParent != null) {
            oldParent.removeChild(chosenConv.getNode2());
            oldParent.addChild(parent);
        }

        if (parent.isRoot())
            acg.setRoot(parent);
        else {
            if (sister.isRoot())
                acg.setRoot(sister);
        }

        logHGF += getAffectedRegionProb(chosenConv) + getEdgeAttachmentProb(chosenConv);

        // Remove conversion
        acg.deleteConversion(chosenConv);

        // Reject if move has given us a conversion departing from the
        // root edge:
        for (Conversion conv : acg.getConversions()) {
            if (conv.getNode1().isRoot())
                return Double.NEGATIVE_INFINITY;
        }

        return logHGF;
    }

    /**
     * Replaces a portion of the clonal frame with a conversion.
     *
     * @return log Hastings-Green factor
     */
    public double createConversion() {
        double logHGF = 0.0;

        Conversion newConv = new Conversion();

        // Choose affected sites:
        logHGF -= drawAffectedRegion(newConv);

        // Choose attchment points:
        logHGF -= attachEdge(newConv);

        // Skip invisible conversions:
        if (newConv.getNode1() == newConv.getNode2()) {
            return Double.NEGATIVE_INFINITY;
        }

        // Check for conversions which attach above chosen point
        for (Conversion conv : acg.getConversions()) {
            if ((conv.getNode1() == newConv.getNode1()
                && conv.getHeight1() > newConv.getHeight1())
                || (conv.getNode2() == newConv.getNode1()
                && conv.getHeight2() > newConv.getHeight1())) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        Node parent = newConv.getNode1().getParent();
        Node sister = getSibling(newConv.getNode1());
        double newHeight2 = parent.getHeight();
        Node newNode2 = sister;

        // Detach parent from original location above sister
        for (Conversion conv : acg.getConversions()) {
            if (conv.getNode1() == parent) {
                conv.setNode1(sister);
            }

            if (conv.getNode2() == parent) {
                conv.setNode2(sister);
            }
        }

        if (newConv.getNode2() == parent) {
            newConv.setNode2(sister);
        }

        parent.removeChild(sister);
        if (parent.isRoot()) {
            sister.setParent(null);
        } else {
            Node grandParent = parent.getParent();
            grandParent.removeChild(parent);
            grandParent.addChild(sister);
            parent.setParent(null);
        }

        // Attach parent to new location above newConv.node2
        parent.setHeight(newConv.getHeight2());

        if (!newConv.getNode2().isRoot()) {
            Node grandParent = newConv.getNode2().getParent();
            grandParent.removeChild(newConv.getNode2());
            grandParent.addChild(parent);
        }
        parent.addChild(newConv.getNode2());

        for (Conversion conv : acg.getConversions()) {
            if ((conv.getNode1() == newConv.getNode2())
                && (conv.getHeight1() > parent.getHeight())) {
                conv.setNode1(parent);
            }

            if ((conv.getNode2() == newConv.getNode2())
                && (conv.getHeight2() > parent.getHeight())) {
                conv.setNode2(parent);
            }
        }

        // Update newConv so that it replaces the original topology.
        newConv.setHeight2(newHeight2);

        if (!newNode2.isRoot() && newNode2.getParent().getHeight() < newHeight2) {
            newConv.setNode2(newNode2.getParent());
        } else {
            newConv.setNode2(newNode2);
        }

        acg.addConversion(newConv);

        // Ensure root is correctly set:
        // (Can only do this once the tree is stitched back together.)

        if (parent.isRoot())
            acg.setRoot(parent);
        else {
            if (sister.isRoot())
                acg.setRoot(sister);
        }

        // Reject if move has given us a conversion departing from the
        // root edge:
        for (Conversion conv : acg.getConversions()) {
            if (conv.getNode1().isRoot())
                return Double.NEGATIVE_INFINITY;
        }

        // Include reverse move HGF contribution:
        logHGF += Math.log(1.0 / acg.getConvCount());

        return logHGF;
    }

    public static void main(String[] args) throws Exception {

        //Randomizer.setSeed(7);
        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealParameter("1.0"));

        SimulatedACG acg = new SimulatedACG();
        acg.initByName(
            "rho", 0.0001,
            "delta", 500.0,
            "populationModel", popFunc,
            "nTaxa", 10,
            "sequenceLength", 10000);

        /*
         ConversionGraph acg = new ConversionGraph();
         acg.initByName(
         "sequenceLength", 10000,
         //            "fromString", "(0:1.0,1:1.0)2:0.0;");
         "fromString", "[&0,500,0.2,1,800,0.8] (0:1.0,1:1.0)2:0.0;");
         */
        ClonalFrameConversionSwap operator = new ClonalFrameConversionSwap();
        operator.initByName(
            "weight", 1.0,
            "acg", acg,
            "populationModel", popFunc,
            "delta", new RealParameter("50.0"));

        double logHR1, logHR2;

        System.out.println(acg.getExtendedNewick(true));
        do {
            logHR1 = operator.createConversion();
        } while (Double.isInfinite(logHR1));

        System.out.println(acg.getExtendedNewick(true));

        do {
            logHR2 = operator.deleteConversion();
        } while (Double.isInfinite(logHR2));

        System.out.println(acg.getExtendedNewick(true));

        System.out.println("logHR1 = " + logHR1);
        System.out.println("logHR2 = " + logHR2);

    }

}
