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

package argbeast.model;

import argbeast.Recombination;
import argbeast.RecombinationGraph;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import java.util.List;
import java.util.Map;
import org.junit.Test;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SimulatedAlignmentTest {
    
    @Test
    public void test() throws Exception {
        
        Randomizer.setSeed(1);
        
        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealParameter("1.0"));
        
        RecombinationGraph arg = new SimulatedRecombinationGraph();
        arg.initByName(
                "rho", 5.0,
                "delta", 1000.0,
                "populationModel", popFunc,
                "nTaxa", 10,
                "sequenceLength", 10000);
        
        System.out.println(arg);

        // Site model:
        JukesCantor jc = new JukesCantor();
        jc.initByName();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName(
                "mutationRate", new RealParameter("0.5"),
                "substModel", jc);

        // Simulate alignment:
        SimulatedAlignment alignment = new SimulatedAlignment();
        alignment.initByName(
                "arg", arg,
                "siteModel", siteModel,
                "outputFileName", "simulated_alignment.nexus",
                "useNexus", true);
        
        // Compute number of segregating sites in each region
        Map<Recombination, Integer> numSS = Maps.newHashMap();
        int numSSrecomb = 0;
        for (Recombination recomb : arg.getRecombinations()) {
            if (recomb == null)
                continue;
            
            numSS.put(recomb, 0);
            for (long i=recomb.getStartLocus(); i<=recomb.getEndLocus(); i++) {
                int leaf0state = alignment.getCounts().get(0).get((int)i);
                for (int leaf=1; leaf<arg.getLeafNodeCount(); leaf++) {
                    if (alignment.getCounts().get(leaf).get((int)i)!=leaf0state) {
                        numSS.put(recomb, numSS.get(recomb)+1);
                        break;
                    }
                }
            }
            numSSrecomb += numSS.get(recomb);
        }
        
        numSS.put(null, 0);
        for (int i=0; i<arg.getSequenceLength(); i++) {
            int leaf0state = alignment.getCounts().get(0).get((int)i);
            for (int leaf=1; leaf<arg.getLeafNodeCount(); leaf++) {
                if (alignment.getCounts().get(leaf).get((int)i)!=leaf0state) {
                    numSS.put(null, numSS.get(null)+1);
                    break;
                }
            }
        }
        numSS.put(null, numSS.get(null)-numSSrecomb);
        
        for (Recombination recomb : arg.getRecombinations())
            System.out.println(arg.getMarginalNewick(recomb));
        
        // Compare UPGMA topologies with true topologies
        // (Should be enough info here for precise agreement)
        for (Recombination recomb : arg.getRecombinations()) {
            Alignment margAlign = createMarginalAlignment(alignment, arg, recomb);
        }
    }
    
    /**
     * Create Alignment object representing alignment of a region
     * corresponding to a single marginal tree.
     * 
     * @param alignment
     * @param arg
     * @param recomb
     * @return
     * @throws Exception 
     */
    private Alignment createMarginalAlignment(Alignment alignment,
            RecombinationGraph arg, Recombination recomb) throws Exception {
        List<Sequence> sequences = Lists.newArrayList();

        for (int leafIdx=0; leafIdx<alignment.getNrTaxa(); leafIdx++) {
            List<Integer> stateSequence;
            
            if (recomb == null) {
                stateSequence =Lists.newArrayList();
                
                int i=0;
                for (int r=0; r<arg.getNRecombs(); r++) {
                    while (i<arg.getRecombinations().get(r+1).getStartLocus()) {
                        stateSequence.add(alignment.getCounts().get(leafIdx).get(i));
                        i += 1;
                    }
                    i=arg.getRecombinations().get(r+1).getEndLocus()+1;
                }


            } else {

                stateSequence = alignment.getCounts().get(leafIdx)
                        .subList(recomb.getStartLocus(), recomb.getEndLocus()+1);
            }
            
            sequences.add(new Sequence(alignment.getTaxaNames().get(leafIdx),
                    alignment.getDataType().state2string(stateSequence)));
        }
        
        Alignment margAlignment = new Alignment(sequences,
                alignment.getDataType().getStateCount(),
                alignment.getDataType().getDescription());
        
        return margAlignment;
    }
}
