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

import bacter.model.unrestricted.ACGCoalescent;
import bacter.model.unrestricted.SimulatedACG;
import beast.core.OperatorSchedule;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.util.Randomizer;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

/**
 * Test for ClonalFrameConversionSwap operator.
 *
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
public class ClonalFrameConversionSwapTest {

    public ClonalFrameConversionSwapTest() {
    }

    @Test
    public void testHR() throws Exception {

        Randomizer.setSeed(13);

        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealParameter("1.0"));

        SimulatedACG acg = new SimulatedACG();
        acg.initByName(
                "rho", 0.0001,
                "delta", 500.0,
                "populationModel", popFunc,
                "nTaxa", 2,
                "sequenceLength", 10000);

        ACGCoalescent target = new ACGCoalescent();
        target.initByName(
                "acg", acg,
                "populationModel", popFunc,
                "rho", new RealParameter("0.0001"),
                "delta", new RealParameter("500.0"));

        State state = new State();
        state.initByName("stateNode", acg);
        state.initialise();
        state.setPosterior(target);

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

        OperatorSchedule operatorSchedule = new OperatorSchedule();
        operatorSchedule.initByName();
        operator.setOperatorSchedule(operatorSchedule);

        double logHR1, logHR2;

        System.out.println(acg.getExtendedNewick(true));
        do {
            logHR1 = operator.createConversion();
        } while (Double.isInfinite(logHR1));

        operator.accept();
        state.store(0);
        state.acceptCalculationNodes();

        System.out.println(acg.getExtendedNewick(true));
        do {
            logHR2 = operator.deleteConversion();
            if (Math.abs(logHR1+logHR2)>1e-10) {
                operator.reject();
                state.restore();
            }
            else
                break;
        } while (true);

        System.out.println(acg.getExtendedNewick(true));
        System.out.println("logHR1 = " + logHR1);
        System.out.println("logHR2 = " + logHR2);

        assertTrue(Math.abs(logHR1 + logHR2) < 1e-10);
    }

}
