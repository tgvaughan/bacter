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
import bacter.operators.ACGOperator;
import bacter.operators.EdgeCreationOperator;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.util.Randomizer;
import java.io.PrintStream;

/**
 * Operator which reversibly deletes a conversion, modifying the CF
 * to match the marginal tree of the original conversion.
 *
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
public class ClonalFrameConversionSwap extends ConversionCreationOperator {

    @Override
    public double proposal() {
        if (Randomizer.nextBoolean())
            return deleteConversion();
        else
            return createConversion();
    }

    public double deleteConversion() {
        double logHGF = 0.0;

        if (acg.getConvCount()==0)
            return Double.NEGATIVE_INFINITY;

        Conversion conv = acg.getConversions().get(
            Randomizer.nextInt(acg.getConvCount()));

        logHGF -= Math.log(1.0/acg.getConvCount());

        // Abort if conversions attach to node1 above height1.
        for (Conversion otherConv : acg.getConversions()) {
            if (otherConv == conv)
                continue;

            if ((otherConv.getNode1() == conv.getNode1()
                && otherConv.getHeight1() > conv.getHeight1())
                || (otherConv.getNode2() == conv.getNode1()
                && otherConv.getHeight2()> conv.getHeight2()))
                return Double.NEGATIVE_INFINITY;
        }

        // Move all conversion attachments from parent to sister.
        Node parent = conv.getNode1().getParent();
        Node sister = getSibling(conv.getNode1());

        for (Conversion otherConv : acg.getConversions()) {
            if (otherConv == conv)
                continue;

            if (otherConv.getNode1() == parent)
                otherConv.setNode1(sister);

            if (otherConv.getNode2() == parent)
                otherConv.setNode2(sister);
        }

        // Detach node1 from parent:
        parent.removeChild(sister);
        if (parent.getParent() != null) {
            Node grandParent = parent.getParent();
            grandParent.removeChild(parent);
            grandParent.addChild(sister);
        } else
            sister.setParent(null);

        // Move conversions from node2 to parent
        parent.setHeight(conv.getHeight2());

        for (Conversion otherConv : acg.getConversions()) {
            if (otherConv == conv)
                continue;

            if (otherConv.getNode1() == conv.getNode2()
                && otherConv.getHeight1() > conv.getHeight2())
                otherConv.setNode1(parent);

            if (otherConv.getNode2() == conv.getNode2()
                && otherConv.getHeight2() > conv.getHeight2())
                otherConv.setNode2(parent);
        }

        // Make final topology changes:
        Node oldParent = conv.getNode2().getParent();
        parent.addChild(conv.getNode2());
        if (oldParent != null) {
            oldParent.removeChild(conv.getNode2());
            oldParent.addChild(parent);
        }

        logHGF += getAffectedRegionProb(conv) + getEdgeAttachmentProb(conv);

        // Remove conversion
        acg.deleteConversion(conv);

        return logHGF;
    }

    public double createConversion() {
        double logHGF = 0.0;

        Conversion newConv = new Conversion();

        // Choose affected sites:
        logHGF -= drawAffectedRegion(newConv);
        
        // Choose attchment points:
        logHGF -= attachEdge(newConv);

        // Check for conversions which attach above chosen point
        for (Conversion conv : acg.getConversions()) {
            if ((conv.getNode1() == newConv.getNode1()
                && conv.getHeight1()>newConv.getHeight1())
                || (conv.getNode2() == newConv.getNode1()
                && conv.getHeight2()>newConv.getHeight1()))
                return Double.NEGATIVE_INFINITY;
        }

        Node parent = newConv.getNode1().getParent();
        Node sister = getSibling(newConv.getNode1());
        double newHeight2 = parent.getHeight();

        for (Conversion conv : acg.getConversions()) {
            if (conv.getNode1() == parent)
                conv.setNode1(sister);

            if (conv.getNode2() == parent)
                conv.setNode2(sister);
        }

        parent.removeChild(sister);
        if (!parent.isRoot()) {
            Node grandParent = parent.getParent();
            grandParent.removeChild(parent);
            grandParent.addChild(sister);
        }

        parent.setHeight(newConv.getHeight2());
        for (Conversion conv : acg.getConversions()) {
            if ((conv.getNode1() == newConv.getNode2())
                && (conv.getHeight1()>parent.getHeight()))
                conv.setNode1(parent);

            if ((conv.getNode2() == newConv.getNode2())
                && (conv.getHeight2()>parent.getHeight()))
                conv.setNode2(parent);
        }

        if (newConv.getNode2().isRoot()) {
            parent.setParent(null);
            acg.setRoot(parent);
            parent.addChild(newConv.getNode2());
        } else {
            Node grandParent = newConv.getNode2().getParent();
            grandParent.removeChild(newConv.getNode2());
            grandParent.addChild(parent);
        }

        newConv.setNode2(sister);
        newConv.setHeight2(newHeight2);

        acg.addConversion(newConv);

        logHGF += Math.log(1.0/acg.getConvCount());

        return logHGF;
    }

    public static void main(String[] args) throws Exception {
        ConversionGraph acg = new ConversionGraph();
        ConstantPopulation popFunc = new ConstantPopulation();


        ClonalFrameConversionSwap operator = new ClonalFrameConversionSwap();
        operator.initByName(
            "weight", 1.0,
            "acg", acg,
            "populationModel", popFunc,
            "delta", new RealParameter("50.0"));
        popFunc.initByName("popSize", new RealParameter("1.0"));

        acg.initByName(
            "sequenceLength", 10000,
            "fromString", "(0:1.0,1:1.0)2:0.0;");

        System.out.println(acg.getExtendedNewick(true));
        operator.createConversion();
        System.out.println(acg.getExtendedNewick(true));
                
    }
    
}
