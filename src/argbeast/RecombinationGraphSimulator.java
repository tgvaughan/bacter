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

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.evolution.tree.Node;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Simulates an ARG - can be used for chain initialization.")
public class RecombinationGraphSimulator extends RecombinationGraph implements StateNodeInitialiser {

    public Input<Double> rhoInput = new Input<Double>("rho",
            "Recombination rate parameter.", Validate.REQUIRED);
    
    public Input<Double> deltaInput = new Input<Double>("delta",
            "Tract length parameter.", Validate.REQUIRED);
    
    public Input<PopulationFunction> popFuncInput = new Input<PopulationFunction>(
            "populationFunction", "Demographic model to use.", Validate.REQUIRED);

    private double rho, delta;
    private PopulationFunction popFunc;
    
    public RecombinationGraphSimulator() { };
    
    @Override
    public void initAndValidate() throws Exception {

        super.initAndValidate();
        
        rho = rhoInput.get();
        delta = deltaInput.get();
        popFunc = popFuncInput.get();
    }

    /**
     * Use coalescent model to simulate clonal frame.
     */
    private void simulateClonalFrame() {

        // Initialize leaf nodes
        List<Node> leafNodes = Lists.newArrayList();
        for (int i=0; i<alignmentInput.get().getNrTaxa(); i++) {
            Node leaf = new Node();
            leaf.setNr(index);
            if (timeTraitSet != null)
                leaf.setHeight(timeTraitSet.getValue(i));
            else
                leaf.setHeight(0.0);
        }
        
        // Create and sort list of inactive nodes
        List<Node> inactiveNodes = Lists.newArrayList(leafNodes);
        Collections.sort(inactiveNodes, new Comparator<Node>() {
            @Override
            public int compare(Node o1, Node o2) {
                if (o1.getHeight()<o2.getHeight())
                    return -1;
                
                if (o2.getHeight()>o1.getHeight())
                    return 1;
                
                return 0;
            }
        });
        
        List<Node> activeNodes = Lists.newArrayList();
        
        double tau = 0.0;
        while (true) {
            
            // Calculate coalescence propensity
            int k = activeNodes.size();
            double chi = 0.5*k*(k-1);
            tau += Randomizer.nextExponential(chi);
            
            double t = popFunc.getInverseIntensity(tau);
            
            if (inactiveNodes.isEmpty() && activeNodes.size()<2)
                break;
        }
    }
    
    /**
     * Uses a discrete two-state Markov process to generate the list of
     * converted segments along the sequence.
     */
    private void generateRecombinations() {
        
        // Choose the starting state:
    }
    
    @Override
    public void initStateNodes() throws Exception { }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(this);
    }
    
}
