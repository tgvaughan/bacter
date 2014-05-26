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

package argbeast.operators;

import argbeast.Recombination;
import argbeast.RecombinationGraph;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Randomizer;
import feast.input.In;
import java.util.List;

/**
 * Abstract class of ARG operators that use the clonal origin model as the
 * basis for adding new recombinant edges to an existing RecombinationGraph.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public abstract class EdgeCreationOperator extends RecombinationGraphOperator {
    
    public Input<PopulationFunction> popFuncInput = new In<PopulationFunction>(
            "populationModel", "Model of population size through time.").setRequired();

    protected PopulationFunction popFunc;
    
    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
        
        popFunc = popFuncInput.get();
    }
    
    /**
     * Attach chosen recombination to the clonal frame.  Note that only the
     * attachment points (nodes and heights) are set, the affected region of
     * the alignment is not modified.
     * 
     * @param recomb
     * @return log probability density of chosen attachment.
     */
    public double attachEdge(Recombination recomb) {
        
        double logP = 0.0;
        
        List<RecombinationGraph.Event> eventList = arg.getCFEvents();
        
        // Select departure point
        double u = Randomizer.nextDouble()*arg.getClonalFrameLength();
        logP += Math.log(1.0/arg.getClonalFrameLength());
        
        for (Node node : arg.getNodesAsArray()) {
            if (node.isRoot())
                continue;
            
            if (u<node.getLength()) {
                recomb.setNode1(node);
                recomb.setHeight1(node.getHeight() + u);
                break;
            } else
                u -= node.getLength();
        }
        
        // Select arrival point
        logP += coalesceEdge(recomb);
        
        return logP;
    }
    
    /**
     * Retrieve probability density for both attachment points of the given
     * recombinant edge.
     * 
     * @param recomb
     * @return log probability density
     */
    public double getEdgeAttachmentProb(Recombination recomb) {
        double logP = 0.0;
        
        logP += Math.log(1.0/arg.getClonalFrameLength());
        logP += getEdgeCoalescenceProb(recomb);
        
        return logP;
    }
    
    /**
     * Take a recombination with an existing departure point and determine
     * the arrival point by allowing it to coalesce with the clonal frame.
     * 
     * @param recomb recombination to modify
     * @return log probability density of coalescent point chosen.
     */
    public double coalesceEdge(Recombination recomb) {
        double logP = 0.0;
        
        List<RecombinationGraph.Event> events = arg.getCFEvents();
        
        // Locate event immediately below departure point
        int startIdx = 0;
        while (events.get(startIdx+1).getHeight()<recomb.getHeight1())
            startIdx += 1;
                
        // Choose edge length in dimensionless time.
        double u = Randomizer.nextExponential(1.0);
        
        // Determine arrival point in real time
        for (int i=startIdx; i<events.size(); i++) {
            
            RecombinationGraph.Event event = events.get(i);
            
            double t = Math.max(recomb.getHeight1(), event.getHeight());
        
            // Determine length of interval in dimensionless time
            double intervalArea;
            if (i<events.size()-1)
                intervalArea = popFunc.getIntegral(t, events.get(i+1).getHeight());
            else
                intervalArea = Double.POSITIVE_INFINITY;
            
            // Check whether arrival falls within this interval
            if (u<intervalArea*event.getLineageCount()) {
                
                // Set arrival point in real time
                recomb.setHeight2(popFunc.getInverseIntensity(
                        popFunc.getIntensity(t) + u/event.getLineageCount()));
                
                // Attach to random clonal frame lineage extant at this time
                int z = Randomizer.nextInt(event.getLineageCount());
                for (Node node : arg.getNodesAsArray()) {
                    if (recomb.getHeight2()>node.getHeight() &&
                            (node.isRoot() || recomb.getHeight2()<node.getParent().getHeight())) {
                        if (z==0) {
                            recomb.setNode2(node);
                            break;
                        } else
                            z -= 1;
                    }
                }

                logP += -u + Math.log(1.0/popFunc.getPopSize(recomb.getHeight2()));
                break;
                
            } else {
                u -= intervalArea*event.getLineageCount();
                logP += -intervalArea*event.getLineageCount();
            }
        }
        
        return logP;
    }
    
    /**
     * Get probability density for the arrival time of the given recombinant
     * edge under ClonalOrigin's coalescent model.
     * 
     * @param recomb
     * @return log probability density
     */
    public double getEdgeCoalescenceProb(Recombination recomb) {
        double logP = 0.0;
        
        List<RecombinationGraph.Event> events = arg.getCFEvents();
        
        // Find event immediately below departure point
        int startIdx = 0;
        while (events.get(startIdx+1).getHeight()<recomb.getHeight1())
            startIdx += 1;
        
        // Compute probability of edge length and arrival
        for (int i=startIdx; i<events.size() && events.get(i).getHeight()<recomb.getHeight2(); i++) {           
            double t1 = Math.max(recomb.getHeight1(), events.get(i).getHeight());
            double t2 = recomb.getHeight2();
            if (i<events.size()-1)
                t2 = Math.min(t2, events.get(i+1).getHeight());
        
            double intervalArea = popFunc.getIntegral(t1, t2);
            logP += -intervalArea*events.get(i).getLineageCount();
        }
        
        // Probability of single coalescence event
        logP += Math.log(1.0/popFunc.getPopSize(recomb.getHeight2()));
        
        return logP;
    }
    
}
