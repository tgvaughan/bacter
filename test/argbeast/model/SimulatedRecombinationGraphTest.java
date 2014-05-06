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

package argbeast.model;

import argbeast.Recombination;
import argbeast.model.SimulatedRecombinationGraph;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.util.Randomizer;
import static org.junit.Assert.*;
import org.junit.Test;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SimulatedRecombinationGraphTest {
    
    public SimulatedRecombinationGraphTest() {
    }

    /**
     * Generate a random nucleotide sequence of specified length.
     * 
     * @param length
     * @return String representation of DNA sequence
     */
    private String getSeq(int length) {
        StringBuilder seq = new StringBuilder();
        String alphabet = "GTCA";
        for (int i=0; i<length; i++)
            seq.append(alphabet.charAt((Randomizer.nextInt(4))));
        
        return seq.toString();
    }
    
    @Test
    public void test() throws Exception {
        
        Randomizer.setSeed(42);
        
        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealParameter("1.0"));

        Alignment alignment = new Alignment();
        int length = 10000;
        alignment.initByName(
                "sequence", new Sequence("taxon1", getSeq(length)),
                "sequence", new Sequence("taxon2", getSeq(length)),
                "datatype", "nucleotide"
        );
        
        SimulatedRecombinationGraph rgs = new SimulatedRecombinationGraph();
        
        int nIter = 100000;
        double meanNRecomb = 0.0;
        double meanCoalTime = 0.0;
        double meanDepHeight = 0.0;
        double meanArrHeight = 0.0;
        int totalRecomb = 0;
        
        for (int i=0; i<nIter; i++) {
            rgs.initByName(
                    "rho", 1.0,
                    "delta", 50.0,
                    "populationModel", popFunc,
                    "alignment", alignment);

            meanNRecomb += rgs.getNRecombs();
            meanCoalTime += rgs.getRoot().getHeight();
            for (Recombination recomb : rgs.getRecombinations()) {
                if (recomb == null)
                    continue;
                
                totalRecomb += 1;
                meanDepHeight += recomb.getHeight1();
                meanArrHeight += recomb.getHeight2();
            }
        }
        
        meanNRecomb /= nIter;
        meanCoalTime /= nIter;
        meanDepHeight /= totalRecomb;
        meanArrHeight /= totalRecomb;
        
        System.out.println("meanNRecomb=" + meanNRecomb + " (truth 1.0)");
        System.out.println("meanCoalTime=" + meanCoalTime + " (truth 1.0)");
        System.out.println("meanDepHeight=" + meanDepHeight + " (truth 1.0)");
        System.out.println("meanArrHeight=" + meanArrHeight + " (truth 1.66)");
        
        assertTrue(Math.abs(meanCoalTime-1.0)<1e-2);
        assertTrue(Math.abs(meanNRecomb-1.0)<1e-2);
        assertTrue(Math.abs(meanDepHeight-1.0)<1e-2);
        assertTrue(Math.abs(meanArrHeight-1.66)<1e-2);
        
    }
}
