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
@Description("Simulates an ARG - can be used for chain initialization or for "
        + "sampler validation.")
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
        
        simulateClonalFrame();
        generateRecombinations();
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
            public int compare(Node n1, Node n2) {
                if (n1.getHeight()<n2.getHeight())
                    return -1;
                
                if (n2.getHeight()>n1.getHeight())
                    return 1;
                
                return 0;
            }
        });
        
        List<Node> activeNodes = Lists.newArrayList();
        
        double tau = 0.0;
        int nextNr = leafNodes.size();
        while (true) {
            
            // Calculate coalescence propensity
            int k = activeNodes.size();
            double chi = 0.5*k*(k-1);
            
            // Draw scaled coalescent time
            if (chi>0.0)
                tau += Randomizer.nextExponential(chi);
            else
                tau = Double.POSITIVE_INFINITY;
            
            // Convert to real time
            double t = popFunc.getInverseIntensity(tau);
            
            // If new time takes us past next sample time, insert that sample
            if (!inactiveNodes.isEmpty() && t>inactiveNodes.get(0).getHeight()) {
                Node nextActive = inactiveNodes.remove(0);
                activeNodes.add(nextActive);
                tau = popFunc.getIntensity(nextActive.getHeight());
                continue;
            }
            
            // Coalesce random pair of active nodes.
            Node node1 = activeNodes.remove(Randomizer.nextInt(k));
            Node node2 = activeNodes.remove(Randomizer.nextInt(k-1));
            
            Node parent = new Node();
            parent.addChild(node1);
            parent.addChild(node2);
            parent.setHeight(t);
            parent.setNr(nextNr++);
            
            activeNodes.add(parent);
            
            if (inactiveNodes.isEmpty() && activeNodes.size()<2)
                break;
        }
        
        // Remaining active node is root
        setRoot(activeNodes.get(0));
    }
    
    /**
     * Uses a discrete two-state Markov process to generate the list of
     * converted segments along the sequence.
     */
    private void generateRecombinations() {
        
        double pRec = 0.5*rho*getClonalFrameLength()/getSequenceLength();
        double pTractEnd = 1.0/delta;
        double p0cf = 1.0/(1.0 + pRec*delta);

        long l = 0;
        while (true) {
            long dl;
            
            if (l==0 && Randomizer.nextDouble()<p0cf)
                dl = 0;
            else
                dl = Randomizer.nextGeometric(pRec)+1;

            if (dl>=getSequenceLength())
                break;
            
            Recombination recomb = new Recombination();
            recomb.setStartLocus(dl);
            
            dl += Randomizer.nextGeometric(pTractEnd);
            recomb.setEndLocus(Math.min(dl,getSequenceLength()-1));
            
            addRecombination(recomb);
        }
    }
    
    /**
     * Associates recombinations with the clonal frame, selecting points of
     * departure and coalescence.
     */
    private void associateRecombinationsWithCF() { }
    
    @Override
    public void initStateNodes() throws Exception { }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(this);
    }
    
}
