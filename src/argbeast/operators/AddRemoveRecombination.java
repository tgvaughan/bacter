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
package argbeast.operators;

import argbeast.Recombination;
import argbeast.RecombinationGraph;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which adds and removes recombinations to/from an ARG.")
public class AddRemoveRecombination extends RecombinationGraphOperator {
    
    public Input<RealParameter> rhoInput = new Input<RealParameter>("rho",
            "Recombination rate parameter.", Validate.REQUIRED);
    
    public Input<RealParameter> deltaInput = new Input<RealParameter>("delta",
            "Tract length parameter.", Validate.REQUIRED);
    
    public Input<PopulationFunction> popFuncInput = new Input<PopulationFunction>(
            "populationModel", "A population size model.", Validate.REQUIRED);
    
    private RecombinationGraph arg;
    private PopulationFunction popFunc;

    private enum Type { COALESCENCE, SAMPLE };
    private class Event {
        Type type;
        double realTime;
        double dimensionlessTime = -1.0;
        int lineages = -1;
        Node node;
        
        public Event(Node node) {
            this.node = node;
            if (node.isLeaf())
                this.type = Type.SAMPLE;
            else
                this.type = Type.COALESCENCE;
            
            this.realTime = node.getHeight();
        }
    }
    
    private List<Event> eventList;
    private Map<Node, Event> eventMap;
    
    public AddRemoveRecombination() { };
    
    @Override
    public void initAndValidate() {
        arg = argInput.get();
        popFunc = popFuncInput.get();
        
        eventList = Lists.newArrayList();
        eventMap = Maps.newHashMap();
    };

    @Override
    public double proposal() {
        double logHGF = 0;

        updateEvents();
        
        if (Randomizer.nextDouble()<0.5) {
            
            // Add
            
            logHGF += Math.log(1.0/(arg.getNRecombs()+1));
            logHGF -= drawNewRecomb();
            
            if (!arg.isValid())
                return Double.NEGATIVE_INFINITY;
            
        } else {
            
            // Remove
            
            // Select recombination to remove:
            Recombination recomb = arg.getRecombinations().get(Randomizer.nextInt(arg.getNRecombs())+1);
            
            // Calculate HGF
            logHGF += getRecombProb(recomb);
            logHGF -= Math.log(1.0/arg.getNRecombs());
            
            // Remove recombination
            arg.deleteRecombination(recomb);
        }
        
        return logHGF;
    }
    
    /**
     * Add new recombination to ARG, returning the probability density of the
     * new edge and converted region location.
     * 
     * @return log of proposal density
     */
    private double drawNewRecomb() {
        double logP = 0;
        
        Recombination newRecomb = new Recombination();
        
        // Select starting point on clonal frame
        double u = Randomizer.nextDouble()*arg.getClonalFrameLength();
        
        for (Node node : arg.getNodesAsArray()) {
            u -= node.getLength();
            
            if (u<0) {
                newRecomb.setNode1(node);
                newRecomb.setHeight1(u+node.getParent().getHeight());
            }
        }
        logP += -Math.log(arg.getClonalFrameLength());
        
        // Find event corresponding to node below starting point
        Event event1 = eventMap.get(newRecomb.getNode1());
        
        // Transform starting time to dimensionless time
        double tau1 = event1.dimensionlessTime
                + popFunc.getIntegral(event1.realTime, newRecomb.getHeight1());
        
        // Draw coalescent time with clonal frame
        u = Randomizer.nextExponential(1.0);
        int event1idx = eventList.indexOf(event1);
        
        double tau2 = 0;
        int event2idx;
        for (event2idx=event1idx; event2idx<eventList.size(); event2idx++) {
            Event thisEvent = eventList.get(event2idx);
            double lastTau = Math.max(tau1, thisEvent.dimensionlessTime);
            if (event2idx < eventList.size()-1) {
                Event nextEvent = eventList.get(event2idx+1);
                
                if ((nextEvent.dimensionlessTime-lastTau)*thisEvent.lineages>u) {
                   tau2 = lastTau + u/thisEvent.lineages;
                   break;
                } else {
                    u -= (nextEvent.dimensionlessTime-lastTau)*thisEvent.lineages;
                    logP += -popFunc.getIntegral(thisEvent.realTime, nextEvent.realTime);
                }
            } else
                tau2 = lastTau + u;
        }
        
        Event event2 = eventList.get(event2idx);
        newRecomb.setNode2(event2.node);
        newRecomb.setHeight2(event2.realTime
                + popFunc.getInverseIntensity(tau2)
                -popFunc.getInverseIntensity(event2.realTime));

        logP += -popFunc.getIntegral(event2.realTime, newRecomb.getHeight2());
        
        // Draw location of converted region.  Currently draws start locus 
        // uniformly from among available unconverted loci and draws the tract
        // length from an exponential distribution.  If the end of the tract
        // exceeds the start of the next region, the proposal is rejected.
        
        int convertableLength = 0;
        for (int ridx=0; ridx<arg.getNRecombs(); ridx++) {
            Recombination recomb = arg.getRecombinations().get(ridx);
            
            if (ridx==0) {
                convertableLength += Math.max(0, recomb.getStartLocus()-1);
            } else {
                Recombination prevRecomb = arg.getRecombinations().get(ridx-1);
                convertableLength += Math.max(0,
                        recomb.getStartLocus()-prevRecomb.getEndLocus()-3);
            }
            
            if (ridx==arg.getNRecombs()-1) {
                int finalLocus = arg.getSequenceLength()-1;
                convertableLength += Math.max(0, finalLocus-recomb.getEndLocus()-1);
            }
        }
        
        int z = Randomizer.nextInt(convertableLength);
        logP += Math.log(1.0/convertableLength);
        
        for (int ridx=0; ridx<arg.getNRecombs(); ridx++) {
            Recombination recomb = arg.getRecombinations().get(ridx);
            
            if (ridx==0) {
                z -= Math.max(0, recomb.getStartLocus()-1);
                
                if (z<0) {
                    newRecomb.setStartLocus(Math.max(0, recomb.getStartLocus()-1)+z);
                    break;
                }
            } else {
                Recombination prevRecomb = arg.getRecombinations().get(ridx-1);
                z -= Math.max(0,
                        recomb.getStartLocus()-prevRecomb.getEndLocus()-3);
                
                if (z<0) {
                    newRecomb.setStartLocus(recomb.getStartLocus()-1+z);
                    break;
                }
            }
            
            if (ridx==arg.getNRecombs()-1) {
                int finalLocus = arg.getSequenceLength()-1;
                z -= Math.max(0, finalLocus-recomb.getEndLocus()-1);
                
                if (z<0) {
                    newRecomb.setStartLocus(finalLocus+1+z);
                }
            }
        }
        
        arg.addRecombination(newRecomb);
        
        int convertedLength = (int)Randomizer
                .nextExponential(1.0/deltaInput.get().getValue());
        logP += -convertedLength/deltaInput.get().getValue();
        
        newRecomb.setEndLocus(newRecomb.getStartLocus()+convertedLength);
        
        return logP;
    }
    
    /**
     * Obtain proposal density for the move which results in the current state
     * by adding the recombination recomb to a state without that recombination.
     * 
     * @param recomb
     * @return log of proposal density
     */
    private double getRecombProb(Recombination recomb) {
        double logP = 0;
        
        // Probability of choosing random point on clonal frame
        logP += Math.log(1.0/arg.getClonalFrameLength());
        
        // Find event corresponding to node below starting point
        Event event1 = eventMap.get(recomb.getNode1());
        
        // Transform starting time to dimensionless time
        double tau1 = event1.dimensionlessTime
                + popFunc.getIntegral(event1.realTime, recomb.getHeight1());
        
        // Find event corresponding to node below finishing point
        Event event2 = eventMap.get(recomb.getNode2());
        
        // Transform finishing time to dimensionless time
        double tau2 = event2.dimensionlessTime
                + popFunc.getIntegral(event2.realTime, recomb.getHeight2());
        
        // Calculate probability of coalescent time with clonal frame
        int event1idx = eventList.indexOf(event1);
        int event2idx = eventList.indexOf(event2);
        
        for (int eidx=event1idx; eidx<=event2idx; eidx++) {
            
            double tauA = Math.max(eventList.get(eidx).dimensionlessTime, tau1);
            double tauB;
            if (eidx<eventList.size()-1) {
                tauB = Math.min(eventList.get(eidx+1).dimensionlessTime, tau2);
            } else
                tauB = tau2;
            logP += (tauB-tauA)*eventList.get(eidx).lineages;
        }

        // Calculate probability of converted region.
        int convertableLength = 0;
        for (int ridx=0; ridx<arg.getNRecombs(); ridx++) {
            Recombination thisRecomb = arg.getRecombinations().get(ridx);
            
            if (ridx==0) {
                convertableLength += Math.max(0, thisRecomb.getStartLocus()-1);
            } else {
                Recombination prevRecomb = arg.getRecombinations().get(ridx-1);
                convertableLength += Math.max(0,
                        thisRecomb.getStartLocus()-prevRecomb.getEndLocus()-3);
            }
            
            if (ridx==arg.getNRecombs()-1) {
                int finalLocus = arg.getSequenceLength()-1;
                convertableLength += Math.max(0, finalLocus-thisRecomb.getEndLocus()-1);
            }
        }
        
        logP += Math.log(1.0/convertableLength)
                + -(recomb.getEndLocus()-recomb.getStartLocus())*deltaInput.get().getValue();
        
        return logP;
    }
    
    /**
     * Assemble sorted list of events on clonal frame and a map from nodes
     * to these events.
     */
    private void updateEvents() {
        eventList.clear();
        eventMap.clear();
        
        // Create event list
        for (Node node : arg.getNodesAsArray()) {
            Event event = new Event(node);
            eventList.add(event);
            eventMap.put(node, event);
        }
        
        // Sort events in increasing order of their times
        Collections.sort(eventList, new Comparator<Event>() {
            @Override
            public int compare(Event o1, Event o2) {
                if (o1.realTime<o2.realTime)
                    return -1;
                
                if (o2.realTime<o1.realTime)
                    return 1;
                
                return 0;
            }
        });
        
        // Compute dimensionless event times and lineage counts
        double tau = eventList.get(0).realTime;
        int k = 0;
        
        for (int eidx = 0; eidx<eventList.size(); eidx++) {
            Event event = eventList.get(eidx);            
            
            // Record current dimensionless time
            event.dimensionlessTime = tau;
            
            // Update and record current lineage count
            if (event.type == Type.SAMPLE)
                k += 1;
            else
                k -= 1;
            event.lineages = k;
            
            // Increment dimensionless time using population model
            if (eidx<eventList.size()-1) {
                Event nextEvent = eventList.get(eidx+1);
                tau += popFunc.getIntegral(event.realTime, nextEvent.realTime);
            }
        }
    }
}
