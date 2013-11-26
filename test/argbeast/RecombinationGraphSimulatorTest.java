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

package argbeast;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.util.Randomizer;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class RecombinationGraphSimulatorTest {
    
    public RecombinationGraphSimulatorTest() {
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
        
        Randomizer.setSeed(53);
        
        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealParameter("1.0"));

        Alignment alignment = new Alignment();
        int length = 10000;
        alignment.initByName(
                "sequence", new Sequence("taxon1", getSeq(length)),
                "sequence", new Sequence("taxon2", getSeq(length)),
                "sequence", new Sequence("taxon3", getSeq(length)),
                "sequence", new Sequence("taxon4", getSeq(length)),
                "sequence", new Sequence("taxon5", getSeq(length)),
                "sequence", new Sequence("taxon6", getSeq(length)),
                "sequence", new Sequence("taxon7", getSeq(length)),
                "sequence", new Sequence("taxon8", getSeq(length)),
                "sequence", new Sequence("taxon9", getSeq(length)),
                "sequence", new Sequence("taxon10", getSeq(length)),
                "datatype", "nucleotide"
        );
        
        RecombinationGraphSimulator rgs = new RecombinationGraphSimulator();
        
        RecombinationGraphStats stats = new RecombinationGraphStats();
        stats.initByName("arg", rgs);
        stats.init(System.out);
        System.out.println();
        
        for (int i=0; i<100; i++) {
            rgs.initByName(
                    "rho", 1.0,
                    "delta", 50.0,
                    "populationFunction", popFunc,
                    "alignment", alignment);
            
            stats.log(-1, System.out);
            System.out.println();
        }
    }
}
