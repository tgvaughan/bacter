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

import bacter.*;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.util.ClusterTree;
import beast.util.Randomizer;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

/**
 * Tests the calculation of the ACG likelihood given the sequence data.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ConversionGraphLikelihoodTest extends TestBase {
    
    public ConversionGraphLikelihoodTest() { }
    
    @Test
    public void testClonalFrameLikelihood() throws Exception {

        Locus locus = new Locus(getAlignment());
        locus.setID("locus");

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
                "acg", acg,
                "siteModel", siteModel);
        
        acg.setEverythingDirty(true);
        
        double logP = argLikelihood.calculateLogP();
        double logPtrue = slowLikelihood(acg, locus, locus.getAlignment(), siteModel);
        //double logPtrue = -6444.862402765536;
        
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
        logPtrue = slowLikelihood(acg, locus, locus.getAlignment(), siteModel);
        //logPtrue = -6445.810702954902;
        
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
        logPtrue = slowLikelihood(acg, locus, locus.getAlignment(), siteModel);
        //logPtrue = -6452.466389537251;
        
        relativeDiff = Math.abs(2.0*(logPtrue-logP)/(logPtrue+logP));
        
        assertTrue(relativeDiff<1e-14);
    }
    
    @Test
    public void testLikelihoodUsingSimulatedData() throws Exception {

        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealParameter("1.0"));

        Locus locus = new Locus(10000);
        locus.setID("locus");
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
                "alignment", alignment,
                "acg", acg,
                "siteModel", siteModel);
        
        double logP = argLikelihood.calculateLogP();

        // Compare product of likelihoods of "marginal alignments" with
        // likelihod computed using RGL.
        double logPprime = slowLikelihood(acg, locus, alignment, siteModel);
        
        double relError = 2.0*Math.abs(logP-logPprime)/Math.abs(logP + logPprime);
        System.out.format("logP=%g\nlogPprime=%g\nrelError=%g\n",
                logP, logPprime, relError);
        assertTrue(relError<1e-13);
    }
    
    /**
     * Calculate ACG likelihood using the product of the likelihoods of the
     * marginal trees.
     * 
     * @param acg
     * @param locus
     * @param siteModel
     * @return
     * @throws Exception 
     */
    private double slowLikelihood(ConversionGraph acg, Locus locus, Alignment alignment,
            SiteModel siteModel) throws Exception {

        double logP = 0.0;
        for (Region region : acg.getRegions(locus)) {
            Alignment margAlign = createMarginalAlignment(alignment, acg, region);
            Tree margTree = new Tree(new MarginalTree(acg, region).getRoot());
            TreeLikelihood treeLikelihood = new TreeLikelihood();
            treeLikelihood.initByName(
                    "data", margAlign,
                    "tree", margTree,
                    "siteModel", siteModel);
            
            logP += treeLikelihood.calculateLogP();
        }
        
        return logP;
    }
}
