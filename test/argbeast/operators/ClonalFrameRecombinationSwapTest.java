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

package argbeast.operators;

import argbeast.Recombination;
import argbeast.RecombinationGraph;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.util.Randomizer;
import beast.util.TreeParser;
import static org.junit.Assert.*;
import org.junit.Test;

/**
 * Tests for AddRemoveRecombination operator.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ClonalFrameRecombinationSwapTest {
    
    public ClonalFrameRecombinationSwapTest() { }
    
    /**
     * Tests that probability density of forward move calculated
     * by drawNewRecomb() matches probability density of backward
     * move calculated by getRecombProb().
     * 
     * @throws Exception 
     */
    @Test
    public void testHR() throws Exception {
        
        // DEBUG:
        Randomizer.setSeed(42);
        
        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealParameter("1.0"));
        
        TreeParser clonalFrame = new TreeParser("((0:0.5,1:0.5)3:0.2,2:0.7)4:0.0");
        RecombinationGraph arg = new RecombinationGraph();
        arg.assignFrom(clonalFrame);
        arg.initByName("sequenceLength", 10000);
        
        State state = new State();
        state.initByName("stateNode", arg);
        state.initialise();

        Recombination recomb = new Recombination();
        recomb.setStartSite(200);
        recomb.setEndSite(500);
        recomb.setNode1(arg.getNode(0));
        recomb.setNode2(arg.getNode(1));
        recomb.setHeight1(0.1);
        recomb.setHeight2(0.4);
        arg.addRecombination(recomb);

        state.store(0);
        
        // Define operator
        ClonalFrameRecombinationSwap operator = new ClonalFrameRecombinationSwap();
        operator.initByName(
                "weight", 1.0,
                "arg", arg,
                "populationModel", popFunc);

        System.out.println(arg);
        
        // DEBUG (so I don't forget this test is incomplete+broken)
        assertTrue(false);
        
        while (true) {
            double logP1 = operator.proposal();
            System.out.println(arg);
            System.out.println("logP1 = " + logP1);
            
            double logP2 = operator.proposal();
            System.out.println(arg);
            System.out.println("logP2 = " + logP2);
            
            if (Double.isInfinite(logP2))
                break; // WTF?
            
        }
    }
}
