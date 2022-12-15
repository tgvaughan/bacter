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

package bacter.devutils;

import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.Locus;
import beast.base.evolution.tree.Node;
import beast.base.inference.MCMC;
import beast.base.inference.State;
import beast.base.parser.XMLParser;

import java.io.File;
import java.io.PrintStream;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class CFConvSwapExperiment {

    public static Node getSibling(Node node) {
        Node parent = node.getParent();
        if (parent.getLeft() == node)
            return parent.getRight();
        else
            return parent.getLeft();
    }

    public static void disconnectEdge(ConversionGraph acg, Node node) {

        if (node.isRoot())
            throw new IllegalArgumentException("Programmer error: "
                    + "root argument passed to disconnectEdge().");

        Node parent = node.getParent();
        Node sister = getSibling(node);

        if (parent.isRoot()) {
            parent.removeChild(sister);
            sister.setParent(null);
        } else {
            Node grandParent = parent.getParent();
            grandParent.removeChild(parent);
            parent.setParent(null);
            parent.removeChild(sister);
            grandParent.addChild(sister);
        }

        for (Locus locus : acg.getConvertibleLoci()) {
            for (Conversion conv : acg.getConversions(locus)) {
                if (conv.getNode1() == parent)
                    conv.setNode1(sister);

                if (conv.getNode2() == parent)
                    conv.setNode2(sister);
            }
        }
    }

    public static void connectEdge(ConversionGraph acg, Node node, Node destEdgeBase, double destTime) {

        if (node.isRoot())
            throw new IllegalArgumentException("Programmer error: "
                    + "root argument passed to connectEdge().");

        Node parent = node.getParent();

        if (destEdgeBase.isRoot()) {
            parent.addChild(destEdgeBase);
        } else {
            Node grandParent = destEdgeBase.getParent();
            grandParent.removeChild(destEdgeBase);
            grandParent.addChild(parent);
            parent.addChild(destEdgeBase);
        }

        parent.setHeight(destTime);

        for (Locus locus : acg.getConvertibleLoci()) {
            for (Conversion conv : acg.getConversions(locus)) {
                if (conv.getNode1() == destEdgeBase && conv.getHeight1() > destTime)
                    conv.setNode1(parent);

                if (conv.getNode2() == destEdgeBase && conv.getHeight2() > destTime)
                    conv.setNode2(parent);
            }
        }
    }

    public static void main(String[] args) throws Exception {


        // Load model
        XMLParser parser = new XMLParser();
        MCMC mcmc = (MCMC)parser.parseFile(new File("inferencePreSimulatedData.xml"));

        State state = mcmc.startStateInput.get();
        state.setStateFileName("problem.state");
        state.restoreFromFile();

        double oldPosterior = state.robustlyCalcPosterior(mcmc.posteriorInput.get());

        ConversionGraph acg = (ConversionGraph)state.getStateNode(0);

        PrintStream ps = new PrintStream("proposal.trees");
        ps.println(acg);

        Node srcNode = acg.getNode(3);
        Node srcNodeS = getSibling(srcNode);
        Node destNode = acg.getNode(6);
        double t_srcNodeP = srcNode.getParent().getHeight();

        disconnectEdge(acg, srcNode);


        Locus locus = acg.getConvertibleLoci().get(0);
        Conversion convToReplace = acg.getConversions(locus).get(27);
        acg.deleteConversion(convToReplace);


        Node srcNodeP = srcNode.getParent();

        connectEdge(acg, srcNode, destNode, convToReplace.getHeight2());

        // Move Conversions
        for (Conversion conv : acg.getConversions(locus)) {
            boolean moved = false;

            if (conv.getNode1()==srcNode && conv.getHeight1()>srcNodeP.getHeight()) {
                conv.setNode1(destNode);
                moved = true;
            }

            if (conv.getNode2()==srcNode && conv.getHeight2()>srcNodeP.getHeight()) {
                conv.setNode2(destNode);
                moved = true;
            }

            if (moved) {
                while (conv.getHeight1()>conv.getNode1().getParent().getHeight())
                    conv.setNode1(conv.getNode1().getParent());

                while (conv.getHeight2()>conv.getNode2().getParent().getHeight())
                    conv.setNode2(conv.getNode2().getParent());
            }

        }

        Conversion convNew = new Conversion();
        convNew.setLocus(locus);
        convNew.setStartSite(0);
        convNew.setEndSite(4000);
        convNew.setNode1(srcNode);
        convNew.setNode2(srcNodeS);
        convNew.setHeight1(convToReplace.getHeight1());
        convNew.setHeight2(t_srcNodeP);
        acg.addConversion(convNew);

        ps.println(acg);
        double newPosterior = state.robustlyCalcPosterior(mcmc.posteriorInput.get());

        System.out.println(newPosterior-oldPosterior);

//        state.setStateFileName("proposal.state");
//        state.storeToFile(1);

        // Open state file


    }
}
