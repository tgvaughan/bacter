/*
 * Copyright (C) 2014 Tim Vaughan <tgvaughan@gmail.com>
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

import bacter.model.SimulatedAlignment;
import bacter.model.SimulatedRecombinationGraph;
import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.TestBase;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.util.ClusterTree;
import beast.util.Randomizer;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SimulatedAlignmentTest extends TestBase {
    
    @Test
    public void test() throws Exception {
        
        Randomizer.setSeed(5);
        
        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealParameter("1.0"));
        
        ConversionGraph arg = new SimulatedRecombinationGraph();
        arg.initByName(
                "rho", 5.0/10000,
                "delta", 1000.0,
                "populationModel", popFunc,
                "nTaxa", 5,
                "sequenceLength", 10000);
        
        System.out.println(arg);

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
                "arg", arg,
                "siteModel", siteModel,
                "outputFileName", "simulated_alignment.nexus",
                "useNexus", true);
        

        // Compare UPGMA topologies with true topologies
        // (Should be enough info here for precise agreement)
        for (Conversion recomb : arg.getConversions()) {
            Alignment margAlign = createMarginalAlignment(alignment, arg, recomb);
            
            ClusterTree upgmaTree = new ClusterTree();
            upgmaTree.initByName(
                    "clusterType", "upgma",
                    "taxa", margAlign);

            assertTrue(topologiesEquivalent(arg.getMarginalTree(recomb, null), upgmaTree));
        }
    }
    

    
}
