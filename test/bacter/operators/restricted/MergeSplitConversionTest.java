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

package bacter.operators.restricted;

import bacter.operators.restricted.MergeSplitConversion;
import bacter.Conversion;
import bacter.ConversionGraph;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.util.TreeParser;
import static org.junit.Assert.*;
import org.junit.Test;

/**
 * Tests for AddRemoveRecombination operator.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class MergeSplitConversionTest {
    
    public MergeSplitConversionTest() { }
    
    /**
     * Tests that probability density of forward move calculated
     * by drawNewRecomb() matches probability density of backward
     * move calculated by getRecombProb().
     * 
     * @throws Exception 
     */
    @Test
    public void testHR() throws Exception {
        
        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealParameter("1.0"));
        
        TreeParser clonalFrame = new TreeParser("((0:0.5,1:0.5):0.2,2:0.7):0.0");
        ConversionGraph acg = new ConversionGraph();
        acg.assignFrom(clonalFrame);
        acg.initByName("sequenceLength", 10000);
        
        State state = new State();
        state.initByName("stateNode", acg);
        state.initialise();

        Conversion recomb = new Conversion();
        recomb.setStartSite(200);
        recomb.setEndSite(500);
        recomb.setNode1(acg.getNode(0));
        recomb.setNode2(acg.getNode(1));
        recomb.setHeight1(0.1);
        recomb.setHeight2(0.4);
        acg.addConversion(recomb);

        state.store(0);
        
        // Define operator
        MergeSplitConversion operator = new MergeSplitConversion();
        operator.initByName("weight", 1.0,
                "acg", acg,
                "expectedGapSize", "10");

        // Loop until a valid proposal is made                    
        double logP1;
        do {
            logP1 = operator.splitProposal();
            
            if (Double.isInfinite(logP1))
                acg.restore();
            
        } while(Double.isInfinite(logP1));
        
        System.out.println("logP1 = " + logP1);
        
        double logP2 = operator.mergeProposal();
        System.out.println("logP2 = " + logP2);
        
        assertTrue(Math.abs(logP1+logP2)<1e-13);
    }
}
