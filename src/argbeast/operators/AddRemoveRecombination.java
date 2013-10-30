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
        
        Recombination recomb = new Recombination();
        
        // Select starting point on clonal frame
        double u = Randomizer.nextDouble()*arg.getClonalFrameLength();
        
        for (Node node : arg.getNodesAsArray()) {
            u -= node.getLength();
            
            if (u<0) {
                recomb.setNode1(node);
                recomb.setHeight1(u+node.getParent().getHeight());
            }
        }
        
        // Find event corresponding to node below starting point
        Event event1 = eventMap.get(recomb.getNode1());
        
        // Transform starting time to dimensionless time
        double tau1 = event1.dimensionlessTime
                + popFunc.getIntegral(event1.realTime, recomb.getHeight1());
        
        // Draw coalescent time with clonal frame
        
        
        // Draw location of converted region
        
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
