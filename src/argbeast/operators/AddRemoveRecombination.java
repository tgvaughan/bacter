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
import argbeast.model.SimulatedRecombinationGraph;
import argbeast.util.RecombinationGraphStatsLogger;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import java.io.PrintStream;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

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
    
    private enum EventType {COAL, SAMP};
    private class Event {
        EventType type;
        double t;
        double tau = -1.0;
        int lineages = -1;
        
        public Event(Node node) {
            if (node.isLeaf())
                type = EventType.SAMP;
            else
                type = EventType.COAL;
            
            t = node.getHeight();
        }
        
        @Override
        public String toString() {
            return String.format("%s (k:%d t:%g tau:%g)",
                    type.name(), lineages, t, tau);
        }
    }
    
    private List<Event> eventList;
    
    public AddRemoveRecombination() { };
    
    @Override
    public void initAndValidate() {
        arg = argInput.get();
        popFunc = popFuncInput.get();
        
        eventList = Lists.newArrayList();
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
            
            if (arg.getNRecombs()==0)
                return Double.NEGATIVE_INFINITY;
            
            // Select recombination to remove:
            Recombination recomb = arg.getRecombinations().get(
                    Randomizer.nextInt(arg.getNRecombs())+1);
            
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
    public double drawNewRecomb() {
        double logP = 0;
        
        Recombination newRecomb = new Recombination();
        
        // Select starting point on clonal frame
        double u = Randomizer.nextDouble()*arg.getClonalFrameLength();
        logP += Math.log(1.0/arg.getClonalFrameLength());
        
        double tauStart = 0.0;
        
        boolean started = false;
        for (int eidx=0; eidx<eventList.size(); eidx++) {
            Event event = eventList.get(eidx);
            
            if (!started) {
                
                // If edge hasn't started, eidx must be < eventList.size()-1
                double interval = eventList.get(eidx+1).t - event.t;
                
                if (u<interval*event.lineages) {
                    for (Node node : arg.getNodesAsArray()){
                        if (node.isRoot())
                            continue;
                        
                        if (node.getHeight()<=event.t
                                && node.getParent().getHeight()>event.t) {
                            if (u<interval) {
                                newRecomb.setNode1(node);
                                newRecomb.setHeight1(event.t + u);
                                tauStart = popFunc.getIntensity(event.t + u);
                                break;
                            } else
                                u -= interval;
                        }
                    }

                    started = true;
                    u = Randomizer.nextExponential(1.0);
                } else
                    u -= interval*event.lineages;
            }
            
            if (started) {
                double interval;
                if (eidx<eventList.size()-1)
                    interval = (eventList.get(eidx+1).tau
                            - Math.max(event.tau, tauStart));
                else
                    interval = Double.POSITIVE_INFINITY;

                if (u < interval*event.lineages) {
                                        
                    // Fix height of attachment point
                    double tauEnd = Math.max(event.tau, tauStart) + u/event.lineages;
                    double tEnd = popFunc.getInverseIntensity(tauEnd);
                    newRecomb.setHeight2(tEnd);
                    logP += -u
                            + Math.log(1.0/popFunc.getPopSize(tEnd));
                    
                    // Choose particular lineage to attach to
                    int nodeNumber = Randomizer.nextInt(event.lineages);
                    for (Node node : arg.getNodesAsArray()) {
                        if (node.getHeight()<=event.t
                                && (node.isRoot() || node.getParent().getHeight()>event.t)) {
                            if (nodeNumber==0) {
                                newRecomb.setNode2(node);
                                break;
                            } else
                                nodeNumber -= 1;
                        }
                    }
                    break;
                } else {
                    u -= interval*event.lineages;
                    logP += -interval*event.lineages;
                }
                
            }
        }
        
        // Draw location of converted region.  Currently draws start locus 
        // uniformly from among available unconverted loci and draws the tract
        // length from an exponential distribution.  If the end of the tract
        // exceeds the start of the next region, the proposal is rejected.
        
        int convertableLength = arg.getSequenceLength();
        for (int ridx=0; ridx<arg.getNRecombs(); ridx++) {
            Recombination recomb = arg.getRecombinations().get(ridx+1);
            
            convertableLength -= recomb.getEndLocus()-recomb.getStartLocus()+3;
            
            if (recomb.getStartLocus()==0)
                convertableLength += 1;
            
            if (recomb.getEndLocus()==arg.getSequenceLength()-1)
                convertableLength += 1;
        }
        
        int z = Randomizer.nextInt(convertableLength);
        logP += Math.log(1.0/convertableLength);
        
        for (int ridx=0; ridx<arg.getNRecombs(); ridx++) {
            Recombination recomb = arg.getRecombinations().get(ridx+1);
            
            if (z<recomb.getStartLocus()-1)
                break;
            
            z += recomb.getEndLocus()-Math.max(0,recomb.getStartLocus()-1);
        }
        
        newRecomb.setStartLocus(z);
        
        int convertedLength = (int)Randomizer
                .nextGeometric(1.0/deltaInput.get().getValue());
        logP += convertedLength*Math.log(1.0-1.0/deltaInput.get().getValue())
                + Math.log(1.0/deltaInput.get().getValue());
                        
        newRecomb.setEndLocus(newRecomb.getStartLocus()+convertedLength);

        if (!arg.addRecombination(newRecomb))
            return Double.NEGATIVE_INFINITY;
        
        return logP;
    }
      
    /**
     * Obtain proposal density for the move which results in the current state
     * by adding the recombination recomb to a state without that recombination.
     * 
     * @param recomb
     * @return log of proposal density
     */
    public double getRecombProb(Recombination recomb) {
        double logP = 0;
        
        // Probability of choosing random point on clonal frame
        logP += Math.log(1.0/arg.getClonalFrameLength());
        
        double tau1 = popFunc.getIntensity(recomb.getHeight1());
        double tau2 = popFunc.getIntensity(recomb.getHeight2());
        
        boolean started = false;
        for (int eidx=0; eidx<eventList.size(); eidx++) {
            
            Event event = eventList.get(eidx);
            
            if (!started) {
                if (eventList.get(eidx+1).t>recomb.getHeight1())
                    started = true;
            }
            
            if (started) {
                double tauStart = Math.max(event.tau,tau1);

                if (eidx==eventList.size()-1 || eventList.get(eidx+1).tau>tau2) {
                    logP += -event.lineages*(tau2-tauStart) +
                            Math.log(1.0/popFunc.getPopSize(recomb.getHeight2()));
                    break;
                } else {
                    logP += -event.lineages*(eventList.get(eidx+1).tau-tauStart);
                }
            }
        }
        
        // Calculate probability of converted region.
        int convertableLength = arg.getSequenceLength();
        for (int ridx=0; ridx<arg.getNRecombs(); ridx++) {
            Recombination thisRecomb = arg.getRecombinations().get(ridx+1);
            if (thisRecomb == recomb)
                continue;
            
            convertableLength -= thisRecomb.getEndLocus()-thisRecomb.getStartLocus()+3;
            
            if (thisRecomb.getStartLocus()==0)
                convertableLength += 1;
            
            if (thisRecomb.getEndLocus()==arg.getSequenceLength()-1)
                convertableLength += 1;
        }
        
        logP += Math.log(1.0/convertableLength)
                + (recomb.getEndLocus()-recomb.getStartLocus())
                *Math.log(1.0-1.0/deltaInput.get().getValue())
                + Math.log(1.0/deltaInput.get().getValue());
        
        return logP;
    }
    
    /**
     * Assemble sorted list of events on clonal frame and a map from nodes
     * to these events.
     */
    public void updateEvents() {
        eventList.clear();
        
        // Create event list
        for (Node node : arg.getNodesAsArray()) {
            Event event = new Event(node);
            eventList.add(event);
        }
        
        // Sort events in increasing order of their times
        Collections.sort(eventList, new Comparator<Event>() {
            @Override
            public int compare(Event o1, Event o2) {
                if (o1.t<o2.t)
                    return -1;
                
                if (o2.t<o1.t)
                    return 1;
                
                return 0;
            }
        });
        
        // Compute dimensionless event times and lineage counts
        double tau = eventList.get(0).t;
        int k = 0;
        
        for (int eidx = 0; eidx<eventList.size(); eidx++) {
            Event event = eventList.get(eidx);            
            
            // Record current dimensionless time
            event.tau = tau;
            
            // Update and record current lineage count
            if (event.type == EventType.SAMP)
                k += 1;
            else
                k -= 1;
            event.lineages = k;
            
            // Increment dimensionless time using population model
            if (eidx<eventList.size()-1) {
                Event nextEvent = eventList.get(eidx+1);
                tau += popFunc.getIntegral(event.t, nextEvent.t);
            }
        }
    }

}
