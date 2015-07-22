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

import bacter.CFEventList;
import bacter.Conversion;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Randomizer;
import java.util.List;

/**
 * Abstract class of ACG operators that use the clonal origin model as the basis
 * for adding new recombinant edges to an existing ConversionGraph.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public abstract class EdgeCreationOperator extends ACGOperator {

    public Input<PopulationFunction> popFuncInput = new Input<>(
            "populationModel", "Model of population size through time.",
            Input.Validate.REQUIRED);

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
     * @param conv conversion
     * @return log probability density of chosen attachment.
     */
    public double attachEdge(Conversion conv) {
        
        double logP = 0.0;
        
        // Select departure point
        double u = Randomizer.nextDouble()*acg.getClonalFrameLength();
        logP += Math.log(1.0/acg.getClonalFrameLength());
        
        for (Node node : acg.getNodesAsArray()) {
            if (node.isRoot())
                continue;
            
            if (u<node.getLength()) {
                conv.setNode1(node);
                conv.setHeight1(node.getHeight() + u);
                break;
            } else
                u -= node.getLength();
        }
        
        // Select arrival point
        logP += coalesceEdge(conv);
        
        return logP;
    }
    
    /**
     * Retrieve probability density for both attachment points of the given
     * recombinant edge.
     * 
     * @param conv conversion
     * @return log probability density
     */
    public double getEdgeAttachmentProb(Conversion conv) {
        double logP = 0.0;
        
        logP += Math.log(1.0/acg.getClonalFrameLength());
        logP += getEdgeCoalescenceProb(conv);
        
        return logP;
    }
    
    /**
     * Take a recombination with an existing departure point and determine
     * the arrival point by allowing it to coalesce with the clonal frame.
     * 
     * @param conv recombination to modify
     * @return log probability density of coalescent point chosen.
     */
    public double coalesceEdge(Conversion conv) {
        double logP = 0.0;
        
        List<CFEventList.Event> events = acg.getCFEvents();
        
        // Locate event immediately below departure point
        int startIdx = 0;
        while (events.get(startIdx+1).getHeight()<conv.getHeight1())
            startIdx += 1;
                
        // Choose edge length in dimensionless time.
        double u = Randomizer.nextExponential(1.0);
        
        // Determine arrival point in real time
        for (int i=startIdx; i<events.size(); i++) {
            
            CFEventList.Event event = events.get(i);
            
            double t = Math.max(conv.getHeight1(), event.getHeight());
        
            // Determine length of interval in dimensionless time
            double intervalArea;
            if (i<events.size()-1)
                intervalArea = popFunc.getIntegral(t, events.get(i+1).getHeight());
            else
                intervalArea = Double.POSITIVE_INFINITY;
            
            // Check whether arrival falls within this interval
            if (u<intervalArea*event.getLineageCount()) {
                
                // Set arrival point in real time
                conv.setHeight2(popFunc.getInverseIntensity(
                        popFunc.getIntensity(t) + u/event.getLineageCount()));
                
                // Attach to random clonal frame lineage extant at this time
                int z = Randomizer.nextInt(event.getLineageCount());
                for (Node node : acg.getNodesAsArray()) {
                    if (conv.getHeight2()>node.getHeight() &&
                            (node.isRoot() || conv.getHeight2()<node.getParent().getHeight())) {
                        if (z==0) {
                            conv.setNode2(node);
                            break;
                        } else
                            z -= 1;
                    }
                }

                logP += -u + Math.log(1.0/popFunc.getPopSize(conv.getHeight2()));
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
     * @param conv conversion
     * @return log probability density
     */
    public double getEdgeCoalescenceProb(Conversion conv) {
        double logP = 0.0;
        
        List<CFEventList.Event> events = acg.getCFEvents();
        
        // Find event immediately below departure point
        int startIdx = 0;
        while (events.get(startIdx+1).getHeight()<conv.getHeight1())
            startIdx += 1;
        
        // Compute probability of edge length and arrival
        for (int i=startIdx; i<events.size() && events.get(i).getHeight()<conv.getHeight2(); i++) {           
            double t1 = Math.max(conv.getHeight1(), events.get(i).getHeight());
            double t2 = conv.getHeight2();
            if (i<events.size()-1)
                t2 = Math.min(t2, events.get(i+1).getHeight());
        
            double intervalArea = popFunc.getIntegral(t1, t2);
            logP += -intervalArea*events.get(i).getLineageCount();
        }
        
        // Probability of single coalescence event
        logP += Math.log(1.0/popFunc.getPopSize(conv.getHeight2()));
        
        return logP;
    }
    
}
