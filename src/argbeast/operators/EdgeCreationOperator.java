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
        
        // Select departure point
        double u = Randomizer.nextDouble()*arg.getClonalFrameLength();
        logP += Math.log(1.0/arg.getClonalFrameLength());
        
        for (Node node : arg.getNodesAsArray()) {
            if (node.isRoot())
                continue;
            
            if (u<node.getLength()) {
                recomb.setNode1(node);
                recomb.setHeight1(node.getHeight()+u);
            } else
                u -= node.getLength();
        }
        
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
        
        double u = Randomizer.nextExponential(1.0);
        
        int startEventIdx = 0;
        while (events.get(startEventIdx+1).getHeight()<recomb.getHeight1())
            startEventIdx += 1;
        
        for (int eidx=startEventIdx; eidx<events.size(); eidx++) {
            
            RecombinationGraph.Event event = events.get(eidx);
            
            double t = Math.max(recomb.getHeight1(), event.getHeight());
        
            double intervalArea;
            if (eidx<events.size()-1)
                intervalArea = popFunc.getIntegral(t, events.get(eidx+1).getHeight());
            else
                intervalArea = Double.POSITIVE_INFINITY;
            
            if (u<intervalArea*event.getLineageCount()) {
                
                recomb.setHeight2(popFunc.getIntensity(t) + u/event.getLineageCount());
                
                int z = Randomizer.nextInt(event.getLineageCount());
                for (Node node : arg.getNodesAsArray()) {
                    if (recomb.getHeight2()>=node.getHeight() &&
                            (node.isRoot() || node.getParent().getHeight()>recomb.getHeight2())) {
                        if (z==0) {
                            recomb.setNode2(node);
                            break;
                        } else
                            z -= 1;
                    }
                }
                
                break;
                
            } else {
                logP += -intervalArea*event.getLineageCount();
                u -= intervalArea*event.getLineageCount();
            }
        }
        
        logP += Math.log(1.0/popFunc.getPopSize(recomb.getHeight2()));
        
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
        
        int startEventIdx = 0;
        while (events.get(startEventIdx+1).getHeight()<recomb.getHeight1())
            startEventIdx += 1;
        
        for (int eidx=startEventIdx; eidx<events.size(); eidx++) {           
            double t1 = Math.max(recomb.getHeight1(), events.get(eidx).getHeight());
            double t2 = recomb.getHeight2();
            if (eidx<events.size()-1)
                t2 = Math.min(t2, events.get(eidx+1).getHeight());
        
            double intervalArea = popFunc.getIntegral(t1, t2);
            logP += -intervalArea*events.get(eidx).getLineageCount();
            
            if (eidx==events.size()-1 || events.get(eidx+1).getHeight()>recomb.getHeight2())
                break;
        }
        
        logP += Math.log(1.0/popFunc.getPopSize(recomb.getHeight2()));
        
        return logP;
    }
    
}
