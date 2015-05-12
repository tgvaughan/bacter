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

import bacter.Conversion;
import bacter.TestBase;
import bacter.model.ACGCoalescent;
import bacter.model.SimulatedACG;
import bacter.operators.AddRemoveConversion;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
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
public class AddRemoveConversionTest extends TestBase {
    
    public AddRemoveConversionTest() { }
    
    /**
     * Tests that probability density of forward move calculated
 by drawNewConversion() matches probability density of backward
 move calculated by getConversionProb().
     * 
     * @throws Exception 
     */
    @Test
    public void testHR() throws Exception {
        
        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealParameter("1.0"));

        Alignment alignment = getAlignment();

        SimulatedACG acg = new SimulatedACG();
        acg.initByName(
                "rho", 1.0/alignment.getSiteCount(),
                "delta", 50.0,
                "alignment", alignment,
                "populationModel", popFunc);
        

        AddRemoveConversion operator = new AddRemoveConversion();
        
        // Loop until a valid proposal is made
        double logP1;
        List<Conversion> oldConversions;
        do {
            operator.initByName(
                    "weight", 1.0,
                    "acg", acg,
                    "delta", new RealParameter("50.0"),
                    "populationModel", popFunc);
            
            oldConversions = Lists.newArrayList(
                    acg.getConversions(alignment));
        
            logP1 = operator.drawNewConversion();
        } while (Double.isInfinite(logP1));
        
        System.out.println("logP1 = " + logP1);
        
        // Identify new recomination
        Conversion newRecomb = null;
        for (Conversion recomb : acg.getConversions(alignment)) {
            if (!oldConversions.contains(recomb))
                newRecomb = recomb;
        }
        assertNotNull(newRecomb);
        
        double logP2 = operator.getConversionProb(newRecomb);
        System.out.println("logP2 = " + logP2);
        
        assertTrue(Math.abs(logP1-logP2)<1e-10);
    }

    /**
     * Tests whether probability of proposing a conversion lines up with
     * conversion probability found in ACGCoalescent.
     * @throws java.lang.Exception
     */
    @Test
    public void testProbability() throws Exception {

        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealParameter("1.0"));

        Alignment alignment = getAlignment();
        
        SimulatedACG acg = new SimulatedACG();
        acg.initByName(
                "rho", 1.0/alignment.getSiteCount(),
                "delta", 50.0,
                "alignment", alignment,
                "populationModel", popFunc);
        
        RealParameter rho = new RealParameter(Double.toString(1.0/alignment.getSiteCount()));
        RealParameter delta = new RealParameter("50.0");

        AddRemoveConversion operator = new AddRemoveConversion();
        operator.initByName(
            "weight", 1.0,
            "acg", acg,
            "delta", delta,
            "populationModel", popFunc);

        ACGCoalescent coal = new ACGCoalescent();
        coal.initByName(
            "acg", acg,
            "populationModel", popFunc,
            "rho", rho,
            "delta", delta);
        
        double logP1 = 0.0;
        double logP2 = 0.0;
        for (Conversion conv : acg.getConversions(alignment)) {
            logP1 += operator.getConversionProb(conv);
            logP2 += coal.calculateConversionLogP(conv);
        }

        System.out.println("logP1 = " + logP1);
        System.out.println("logP2 = " + logP2);

        assertTrue(Math.abs(logP1-logP2)<1e-10);

    }
}
