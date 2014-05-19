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
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import feast.input.In;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which adds and removes recombinations to/from an ARG.")
public class AddRemoveRecombination extends RecombinationGraphOperator {
    
    public Input<RealParameter> rhoInput = new In<RealParameter>("rho",
            "Recombination rate parameter.").setRequired();
    
    public Input<RealParameter> deltaInput = new In<RealParameter>("delta",
            "Tract length parameter.").setRequired();
    
    public Input<PopulationFunction> popFuncInput = new In<PopulationFunction>(
            "populationModel", "A population size model.").setRequired();
    
    public Input<TreeIntervals> treeIntervalsInput = new In<TreeIntervals>(
            "treeIntervals", "Tree intervals calculation node.").setRequired();
    
    private PopulationFunction popFunc;
    
    private class ProposalFailed extends Exception { };
    
    public AddRemoveRecombination() { };
    
    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        popFunc = popFuncInput.get();
    };

    @Override
    public double proposal() {
        double logHGF = 0;
        
        if (Randomizer.nextBoolean()) {
            
            // Add
            
            logHGF += Math.log(1.0/(arg.getNRecombs()+1));
            
            try {
                logHGF -= drawNewRecomb();
            } catch (ProposalFailed ex) {
                return Double.NEGATIVE_INFINITY;
            }
            
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
     * @throws argbeast.operators.AddRemoveRecombination.ProposalFailed 
     */
    public double drawNewRecomb() throws ProposalFailed {
        double logP = 0;

        Recombination newRecomb = new Recombination();
        
        List<RecombinationGraph.Event> eventList = arg.getCFEvents();
        
        // Select starting point on clonal frame
        double u = Randomizer.nextDouble()*arg.getClonalFrameLength();
        logP += Math.log(1.0/arg.getClonalFrameLength());
        
        boolean started = false;
        for (int eidx=0; eidx<eventList.size(); eidx++) {
            RecombinationGraph.Event event = eventList.get(eidx);
            
            if (!started) {
                
                // If edge hasn't started, eidx must be < eventList.size()-1
                double interval = eventList.get(eidx+1).getHeight() - event.getHeight();
                
                if (u<interval*event.getLineageCount()) {
                    for (Node node : arg.getNodesAsArray()){
                        if (node.isRoot())
                            continue;
                        
                        if (node.getHeight()<=event.getHeight()
                                && node.getParent().getHeight()>event.getHeight()) {
                            if (u<interval) {
                                newRecomb.setNode1(node);
                                newRecomb.setHeight1(event.getHeight() + u);
                                break;
                            } else
                                u -= interval;
                        }
                    }

                    started = true;
                    u = Randomizer.nextExponential(1.0);
                } else
                    u -= interval*event.getLineageCount();
            }
            
            if (started) {
                double t = Math.max(event.getHeight(), newRecomb.getHeight1());
                double interval;
                if (eidx<eventList.size()-1)
                    interval = popFunc.getIntegral(t, eventList.get(eidx+1).getHeight());
                else
                    interval = Double.POSITIVE_INFINITY;

                if (u < interval*event.getLineageCount()) {
                                        
                    // Fix height of attachment point
                    double tauEnd = popFunc.getIntensity(t) + u/event.getLineageCount();
                    double tEnd = popFunc.getInverseIntensity(tauEnd);
                    newRecomb.setHeight2(tEnd);
                    logP += -u + Math.log(1.0/popFunc.getPopSize(tEnd));
                    
                    // Choose particular lineage to attach to
                    int nodeNumber = Randomizer.nextInt(event.getLineageCount());
                    for (Node node : arg.getNodesAsArray()) {
                        if (node.getHeight()<=event.getHeight()
                                && (node.isRoot() || node.getParent().getHeight()>event.getHeight())) {
                            if (nodeNumber==0) {
                                newRecomb.setNode2(node);
                                break;
                            } else
                                nodeNumber -= 1;
                        }
                    }
                    break;
                } else {
                    u -= interval*event.getLineageCount();
                    logP += -interval*event.getLineageCount();
                }
                
            }
        }
        
        // Draw location of converted region.  Currently draws start locus 
        // uniformly from among available unconverted loci and draws the tract
        // length from an exponential distribution.  If the end of the tract
        // exceeds the start of the next region, the proposal is rejected.
        
        int convertableLength = getConvertibleSiteCount(null);
        if (convertableLength==0)
            throw new ProposalFailed();
        
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
            throw new ProposalFailed();
        
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
        
        List<RecombinationGraph.Event> eventList = arg.getCFEvents();
        
        // Probability of choosing random point on clonal frame
        logP += Math.log(1.0/arg.getClonalFrameLength());
        
        boolean started = false;
        for (int eidx=0; eidx<eventList.size(); eidx++) {
            
            RecombinationGraph.Event event = eventList.get(eidx);
            
            if (!started) {
                if (eventList.get(eidx+1).getHeight()>recomb.getHeight1())
                    started = true;
            }
            
            if (started) {
                double t = Math.max(recomb.getHeight1(), event.getHeight());

                if (eidx==eventList.size()-1
                        || eventList.get(eidx+1).getHeight()>recomb.getHeight2()) {
                    logP += -event.getLineageCount()*popFunc.getIntegral(t, recomb.getHeight2())
                            + Math.log(1.0/popFunc.getPopSize(recomb.getHeight2()));
                    break;
                } else {
                    logP += -event.getLineageCount()
                            *popFunc.getIntegral(t, eventList.get(eidx+1).getHeight());
                }
            }
        }
        
        // Calculate probability of converted region.
        int convertableLength = getConvertibleSiteCount(recomb);
        if (convertableLength==0)
            return Double.NEGATIVE_INFINITY;
        
        logP += Math.log(1.0/convertableLength)
                + (recomb.getEndLocus()-recomb.getStartLocus())
                *Math.log(1.0-1.0/deltaInput.get().getValue())
                + Math.log(1.0/deltaInput.get().getValue());
        
        return logP;
    }
    
    /**
     * Get total number of sites available to convert, optionally skipping
     * one recombination.
     * 
     * @param skip recombination to skip (may be null)
     * @return convertible site count
     */
    private int getConvertibleSiteCount(Recombination skip) {
        int count = 0;
        
        int l=0;
        for (int ridx=1; ridx<=arg.getNRecombs(); ridx++) {
            Recombination recomb = arg.getRecombinations().get(ridx);
            if (recomb == skip)
                continue;
            
            count += Math.max(0, recomb.getStartLocus()-l+1);
            
            l = recomb.getEndLocus()+2;
        }
        
        count += Math.max(0, arg.getSequenceLength() - l); // L - 1 - l + 1
        
        return count;
    }
}
