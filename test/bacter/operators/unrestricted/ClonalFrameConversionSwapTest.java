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

import bacter.model.unrestricted.SimulatedACG;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.ConstantPopulation;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

/**
 * Test for ClonalFrameConversionSwap operator.
 *
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
public class ClonalFrameConversionSwapTest {

    public ClonalFrameConversionSwapTest() { }

    @Test
    public void testHR() throws Exception {

        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealParameter("1.0"));

        SimulatedACG acg = new SimulatedACG();
        acg.initByName(
            "rho", 0.0001,
            "delta", 500.0,
            "populationModel", popFunc,
            "nTaxa", 2,
            "sequenceLength", 10000);

        /*
         ConversionGraph acg = new ConversionGraph();
         acg.initByName(
         "sequenceLength", 10000,
         //            "fromString", "(0:1.0,1:1.0)2:0.0;");
         "fromString", "[&0,500,0.2,1,800,0.8] (0:1.0,1:1.0)2:0.0;");
         */
        ClonalFrameConversionSwap operator = new ClonalFrameConversionSwap();
        operator.initByName(
            "weight", 1.0,
            "acg", acg,
            "populationModel", popFunc,
            "delta", new RealParameter("50.0"));

        double logHR1, logHR2;

        System.out.println(acg.getExtendedNewick(true));
        do {
            logHR1 = operator.createConversion();
        } while (Double.isInfinite(logHR1));

        System.out.println(acg.getExtendedNewick(true));

        do {
            logHR2 = operator.deleteConversion();
        } while (Double.isInfinite(logHR2));

        System.out.println(acg.getExtendedNewick(true));

        System.out.println("logHR1 = " + logHR1);
        System.out.println("logHR2 = " + logHR2);

        assertTrue(Math.abs(logHR1+logHR2)<1e-10);
    }
    
}
