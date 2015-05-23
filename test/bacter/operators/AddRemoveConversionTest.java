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
import bacter.Locus;
import bacter.TestBase;
import bacter.model.ACGCoalescent;
import bacter.model.SimulatedACG;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.coalescent.ConstantPopulation;
import com.google.common.collect.Lists;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

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

        Locus locus = new Locus(10000);
        locus.setID("locus");
        TaxonSet taxonSet = getTaxonSet(10);

        SimulatedACG acg = new SimulatedACG();
        acg.initByName(
                "rho", 1.0/locus.getSiteCount(),
                "delta", 50.0,
                "locus", locus,
                "taxonset", taxonSet,
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
                    acg.getConversions(locus));
        
            logP1 = operator.drawNewConversion();
        } while (Double.isInfinite(logP1));
        
        System.out.println("logP1 = " + logP1);
        
        // Identify new recomination
        Conversion newRecomb = null;
        for (Conversion recomb : acg.getConversions(locus)) {
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

        Locus locus = new Locus(10000);
        locus.setID("locus");

        TaxonSet taxonSet = getTaxonSet(10);
        
        SimulatedACG acg = new SimulatedACG();
        acg.initByName(
                "rho", 1.0/locus.getSiteCount(),
                "delta", 50.0,
                "locus", locus,
                "taxonset", taxonSet,
                "populationModel", popFunc);
        
        RealParameter rho = new RealParameter(Double.toString(1.0/locus.getSiteCount()));
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
        for (Conversion conv : acg.getConversions(locus)) {
            logP1 += operator.getConversionProb(conv);
            logP2 += coal.calculateConversionLogP(conv);
        }

        System.out.println("logP1 = " + logP1);
        System.out.println("logP2 = " + logP2);

        assertTrue(Math.abs(logP1-logP2)<1e-10);

    }
}
