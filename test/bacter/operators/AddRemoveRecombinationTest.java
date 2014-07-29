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

package bacter.operators;

import bacter.operators.AddRemoveRecombination;
import bacter.Recombination;
import bacter.model.SimulatedRecombinationGraph;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.ConstantPopulation;
import com.google.common.collect.Lists;
import java.util.List;
import static org.junit.Assert.*;
import org.junit.Test;

/**
 * Tests for AddRemoveRecombination operator.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class AddRemoveRecombinationTest {
    
    public AddRemoveRecombinationTest() { }
    
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
        
        SimulatedRecombinationGraph arg = new SimulatedRecombinationGraph();
        arg.initByName(
                "rho", 1.0/10000,
                "delta", 50.0,
                "sequenceLength", 10000,
                "nTaxa", 10,
                "populationModel", popFunc);
        

        AddRemoveRecombination operator = new AddRemoveRecombination();
        
        // Loop until a valid proposal is made
        double logP1;
        List<Recombination> oldRecombs;
        do {
            operator.initByName(
                    "weight", 1.0,
                    "arg", arg,
                    "rho", new RealParameter(Double.toString(1.0/10000)),
                    "delta", new RealParameter("50.0"),
                    "populationModel", popFunc);
            
            oldRecombs = Lists.newArrayList(
                    arg.getRecombinations());
        
            logP1 = operator.drawNewRecomb();
        } while (Double.isInfinite(logP1));
        
        System.out.println("logP1 = " + logP1);
        
        // Identify new recomination
        Recombination newRecomb = null;
        for (Recombination recomb : arg.getRecombinations()) {
            if (!oldRecombs.contains(recomb))
                newRecomb = recomb;
        }
        assertNotNull(newRecomb);
        
        double logP2 = operator.getRecombProb(newRecomb);
        System.out.println("logP2 = " + logP2);
        
        assertTrue(Math.abs(logP1-logP2)<1e-10);
    }
}
