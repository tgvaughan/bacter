/*
 * Copyright (C) 2013 Tim Vaughan <tgvaughan@gmail.com>
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
package bacter.model;

import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.Locus;
import bacter.TestBase;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.tree.Node;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.util.ClusterTree;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

/**
 * Tests the calculation of the ACG likelihood given the sequence data.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ACGLikelihoodTest extends TestBase {
    
    public ACGLikelihoodTest() { }
    
    @Test
    public void testClonalFrameLikelihood() throws Exception {

        Locus locus = new Locus("locus", getAlignment());

        // ConversionGraph
        ConversionGraph acg = new ConversionGraph();
        ClusterTree tree = new ClusterTree();
        tree.initByName(
                "clusterType", "upgma",
                "taxa", locus.getAlignment());
        
        acg.assignFrom(tree);
        acg.initByName("locus", locus);
        
        // Site model:
        JukesCantor jc = new JukesCantor();
        jc.initByName();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName(
                "substModel", jc);

        // Likelihood

        ACGLikelihood argLikelihood = new ACGLikelihood();
        argLikelihood.initByName(
                "locus", locus,
                "tree", acg,
                "siteModel", siteModel);


        ACGLikelihoodSlow argLikelihoodSlow = new ACGLikelihoodSlow();
        argLikelihoodSlow.initByName(
                "locus", locus,
                "tree", acg,
                "siteModel", siteModel);

        acg.setEverythingDirty(true);

        double logP = argLikelihood.calculateLogP();
        double logPtrue = argLikelihoodSlow.calculateLogP();

        double relativeDiff = Math.abs(2.0*(logPtrue-logP)/(logPtrue+logP));

        assertTrue(relativeDiff<1e-14);
        
        //Add a single recombination event
        Node node1 = acg.getExternalNodes().get(0);
        Node node2 = node1.getParent();
        double height1 = 0.5*(node1.getHeight() + node1.getParent().getHeight());
        double height2 = 0.5*(node2.getHeight() + node2.getParent().getHeight());
        int startLocus = 100;
        int endLocus = 200;
        Conversion recomb1 = new Conversion(node1, height1, node2, height2,
                startLocus, endLocus, acg, locus);
        acg.addConversion(recomb1);
        
        logP = argLikelihood.calculateLogP();
        logPtrue = argLikelihoodSlow.calculateLogP();

        relativeDiff = Math.abs(2.0*(logPtrue-logP)/(logPtrue+logP));
        
        assertTrue(relativeDiff<1e-14);
        
        // Add another recombination event
        node1 = acg.getExternalNodes().get(0);
        node2 = acg.getNode(20);
        height1 = 0.75*(node1.getHeight() + node1.getParent().getHeight());
        height2 = 0.5*(node2.getHeight() + node2.getParent().getHeight());
        startLocus = 250;
        endLocus = 300;
        Conversion recomb2 = new Conversion(node1, height1, node2, height2,
                startLocus, endLocus, acg, locus);
        acg.addConversion(recomb2);
        
        logP = argLikelihood.calculateLogP();
        logPtrue = argLikelihoodSlow.calculateLogP();

        relativeDiff = Math.abs(2.0*(logPtrue-logP)/(logPtrue+logP));
        
        assertTrue(relativeDiff<1e-14);
    }
    
    @Test
    public void testLikelihoodUsingSimulatedData() throws Exception {

        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealParameter("1.0"));

        Locus locus = new Locus("locus", 10000);
        TaxonSet taxonSet = getTaxonSet(10);
        
        ConversionGraph acg = new SimulatedACG();
        acg.initByName(
                "rho", 5.0/locus.getSiteCount(),
                "delta", 1000.0,
                "populationModel", popFunc,
                "locus", locus,
                "taxonset", taxonSet);
        
        System.out.println(acg);

        // Site model:
        JukesCantor jc = new JukesCantor();
        jc.initByName();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName(
                "mutationRate", new RealParameter("1"),
                "substModel", jc);

        // Simulate alignment:
        SimulatedAlignment alignment = new SimulatedAlignment();
        alignment.initByName(
                "acg", acg,
                "siteModel", siteModel,
                "outputFileName", "simulated_alignment.nexus",
                "useNexus", true);
        
        // Calculate likelihood:
        ACGLikelihood argLikelihood = new ACGLikelihood();
        argLikelihood.initByName(
                "locus", locus,
                "data", alignment,
                "tree", acg,
                "siteModel", siteModel);
        
        double logP = argLikelihood.calculateLogP();

        // Compare product of likelihoods of "marginal alignments" with
        // likelihood computed using RGL.
        ACGLikelihoodSlow argLikelihoodSlow = new ACGLikelihoodSlow();
        argLikelihoodSlow.initByName(
                "locus", locus,
                "data", alignment,
                "tree", acg,
                "siteModel", siteModel);

        double logPprime = argLikelihoodSlow.calculateLogP();


        double relError = 2.0*Math.abs(logP-logPprime)/Math.abs(logP + logPprime);
        System.out.format("logP=%g\nlogPprime=%g\nrelError=%g\n",
                logP, logPprime, relError);
        assertTrue(relError<1e-13);
    }

    @Test
    public void testLikelihoodCaching() throws Exception {
        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealParameter("1.0"));

        Locus locus = new Locus("locus", 10000);
        TaxonSet taxonSet = getTaxonSet(10);


        ConversionGraph acg = new SimulatedACG();
        acg.initByName(
                "rho", 5.0/locus.getSiteCount(),
                "delta", 1000.0,
                "populationModel", popFunc,
                "locus", locus,
                "taxonset", taxonSet);

        State state = new State();
        state.initByName("stateNode", acg);
        state.initialise();

        System.out.println(acg);

        // Site model:
        JukesCantor jc = new JukesCantor();
        jc.initByName();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName(
                "mutationRate", new RealParameter("1"),
                "substModel", jc);

        // Simulate alignment:
        SimulatedAlignment alignment = new SimulatedAlignment();
        alignment.initByName(
                "acg", acg,
                "siteModel", siteModel,
                "outputFileName", "simulated_alignment.nexus",
                "useNexus", true);

        // Calculate likelihood 1:
        ACGLikelihood argLikelihood = new ACGLikelihood();
        argLikelihood.initByName(
                "locus", locus,
                "data", alignment,
                "tree", acg,
                "siteModel", siteModel);

        ACGLikelihoodSlow argLikelihoodSlow = new ACGLikelihoodSlow();
        argLikelihoodSlow.initByName(
                "locus", locus,
                "data", alignment,
                "tree", acg,
                "siteModel", siteModel);

        double logP1 = argLikelihood.calculateLogP();
        double logP1prime = argLikelihoodSlow.calculateLogP();

        double relError = 2.0*Math.abs(logP1-logP1prime)/Math.abs(logP1 + logP1prime);
        System.out.format("logP=%g\nlogPprime=%g\nrelError=%g\n",
                logP1, logP1prime, relError);
        assertTrue(relError<1e-13);

        //Add a single recombination event
        Node node1 = acg.getExternalNodes().get(0);
        Node node2 = node1.getParent();
        double height1 = 0.5*(node1.getHeight() + node1.getParent().getHeight());
        double height2 = 0.5*(node2.getHeight() + node2.getParent().getHeight());
        int startLocus = 500;
        int endLocus = 600;
        Conversion recomb1 = new Conversion(node1, height1, node2, height2,
                startLocus, endLocus, acg, locus);
        acg.addConversion(recomb1);

        double logP2 = argLikelihood.calculateLogP();
        double logP2prime = argLikelihoodSlow.calculateLogP();

        relError = 2.0*Math.abs(logP2-logP2prime)/Math.abs(logP2 + logP2prime);
        System.out.format("logP=%g\nlogPprime=%g\nrelError=%g\n",
                logP2, logP2prime, relError);
        assertTrue(relError<1e-13);

        state.restore();

        double logP3 = argLikelihood.calculateLogP();
        double logP3prime = argLikelihoodSlow.calculateLogP();

        relError = 2.0*Math.abs(logP3-logP3prime)/Math.abs(logP3 + logP3prime);
        System.out.format("logP=%g\nlogPprime=%g\nrelError=%g\n",
                logP3, logP3prime, relError);
        assertTrue(relError<1e-13);
    }

    @Test
    public void testBeagleLikelihood() throws Exception {

        Locus locus = new Locus("locus", getAlignment());

        // ConversionGraph
        ConversionGraph acg = new ConversionGraph();
        ClusterTree tree = new ClusterTree();
        tree.initByName(
                "clusterType", "upgma",
                "taxa", locus.getAlignment());

        acg.assignFrom(tree);
        acg.initByName("locus", locus);

        // Site model:
        JukesCantor jc = new JukesCantor();
        jc.initByName();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName(
                "substModel", jc);

        // Likelihood

        ACGLikelihoodBeagle argLikelihood = new ACGLikelihoodBeagle();
        argLikelihood.initByName(
                "locus", locus,
                "tree", acg,
                "siteModel", siteModel);


        ACGLikelihoodSlow argLikelihoodSlow = new ACGLikelihoodSlow();
        argLikelihoodSlow.initByName(
                "locus", locus,
                "tree", acg,
                "siteModel", siteModel);

        acg.setEverythingDirty(true);

        double logP = argLikelihood.calculateLogP();
        double logPtrue = argLikelihoodSlow.calculateLogP();

        double relativeDiff = Math.abs(2.0*(logPtrue-logP)/(logPtrue+logP));

        assertTrue(relativeDiff<1e-14);

        //Add a single recombination event
        Node node1 = acg.getExternalNodes().get(0);
        Node node2 = node1.getParent();
        double height1 = 0.5*(node1.getHeight() + node1.getParent().getHeight());
        double height2 = 0.5*(node2.getHeight() + node2.getParent().getHeight());
        int startLocus = 100;
        int endLocus = 200;
        Conversion recomb1 = new Conversion(node1, height1, node2, height2,
                startLocus, endLocus, acg, locus);
        acg.addConversion(recomb1);

        logP = argLikelihood.calculateLogP();
        logPtrue = argLikelihoodSlow.calculateLogP();

        relativeDiff = Math.abs(2.0*(logPtrue-logP)/(logPtrue+logP));

        assertTrue(relativeDiff<1e-14);

        // Add another recombination event
        node1 = acg.getExternalNodes().get(0);
        node2 = acg.getNode(20);
        height1 = 0.75*(node1.getHeight() + node1.getParent().getHeight());
        height2 = 0.5*(node2.getHeight() + node2.getParent().getHeight());
        startLocus = 250;
        endLocus = 300;
        Conversion recomb2 = new Conversion(node1, height1, node2, height2,
                startLocus, endLocus, acg, locus);
        acg.addConversion(recomb2);

        logP = argLikelihood.calculateLogP();
        logPtrue = argLikelihoodSlow.calculateLogP();

        relativeDiff = Math.abs(2.0*(logPtrue-logP)/(logPtrue+logP));

        assertTrue(relativeDiff<1e-14);
    }

}
