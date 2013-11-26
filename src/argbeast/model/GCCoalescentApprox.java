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
package argbeast.model;

import argbeast.Recombination;
import argbeast.RecombinationGraph;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.Coalescent;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.evolution.tree.coalescent.TreeIntervals;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Appoximation to the coalescent with gene conversion.")
public class GCCoalescentApprox extends Coalescent {
    
    public Input<RealParameter> rhoInput = new Input<RealParameter>("rho",
            "Recombination rate parameter.", Validate.REQUIRED);
    
    public Input<RealParameter> deltaInput = new Input<RealParameter>("delta",
            "Tract length parameter.", Validate.REQUIRED);
    
    public Input<Boolean> allowSECInput = new Input<Boolean>("allowSameEdgeCoalescence",
    "Allow recombinant edges to attach to the edge they depart from.", true);

    
    RecombinationGraph arg;
    TreeIntervals intervals;
    PopulationFunction popSize;
    int sequenceLength;
    
    boolean allowSEC;
    
    //double logPcfCoal, logPrecombCoal, logPrecombRegion;
    
    @Override
    public void initAndValidate() throws Exception {
        
        intervals = treeIntervalsInput.get();
        if (intervals == null)
            throw new Exception("treeIntervals must be specified.");

        if (!(intervals.treeInput.get() instanceof RecombinationGraph))
            throw new Exception("Error: GCCoalescentApprox only applicable to"
                    + " RecombinationGraph objects.");
        
        arg = (RecombinationGraph)intervals.treeInput.get();
        sequenceLength = arg.getSequenceLength();

        popSize = popSizeInput.get();
        
        allowSEC = allowSECInput.get();
        
        super.initAndValidate();
    }
    
    @Override
    public double calculateLogP() throws Exception {
        
        logP = 0.0;
        
        for (Recombination recomb : arg.getRecombinations())
            logP += calculateRecombinantLogP(recomb);
        
        logP += calculateConvertedRegionMapLogP();
        
        return logP;        
    }

    
    /**
     * Compute probability of recombinant edges under conditional coalescent.
     * @param recomb
     * @return log(P)
     */
    public double calculateRecombinantLogP(Recombination recomb) {

        if (recomb == null) // Clonal frame
            return calculateLogLikelihood(intervals, popSize);
        
        // Probability density of location of recombinant edge start
        double thisLogP = -Math.log(arg.getClonalFrameLength());

        // Identify interval containing the start of the recombinant edge
        int startInterval = 0;
        double t = 0.0;
        while (t+intervals.getInterval(startInterval)<recomb.getHeight1()) {
            t += intervals.getInterval(startInterval);
            startInterval += 1;
        }
        
        // Lineages with which recombinant lineage can coalesce before
        // this time = k(t)-1, while after this time = k(t).
        double oldCoalescenceTime = recomb.getNode1().getParent().getHeight();
        
        int i=startInterval;
        while (t<recomb.getHeight2()) {
            
            double timeA = Math.max(t, recomb.getHeight1());
            
            double timeB;
            int k;
            if (i<intervals.getIntervalCount()) {
                timeB = Math.min(recomb.getHeight2(), t+intervals.getInterval(i));
                k = intervals.getLineageCount(i);
                
                if (!allowSEC && t < oldCoalescenceTime)
                    k -=1;
            } else {
                timeB = recomb.getHeight2();
                k = 1;
            }
            
            double intervalArea = popSize.getIntegral(timeA, timeB);
            thisLogP += -k*intervalArea;
            
            t = timeB;
            i += 1;
        }
        
        // Probability of single coalescence event
        thisLogP += -Math.log(popSize.getPopSize(recomb.getHeight2()));
        
        return thisLogP;
    }
    
    /**
     * Compute probability of number and genome extent of converted segments.
     * @return log(P)
     */
    public double calculateConvertedRegionMapLogP() {
        
        double thisLogP = 0.0;
        
        List<Recombination> recombs = arg.getRecombinations();
        
        double rho = rhoInput.get().getValue();
        double deltaInv = 1.0/deltaInput.get().getValue();
        double cfLength = arg.getClonalFrameLength();
        
        // Probability of recombination per site along sequence
        double alpha = 0.5*rho*cfLength/sequenceLength;
        
        // Probability that sequence begins in the clonal frame:
        double pStartCF = 1.0/(alpha/deltaInv + 1.0);
        
        if (recombs.size()<2){
            // Probability of no recombinations
            thisLogP += Math.log(pStartCF) 
                    + (sequenceLength-1)*Math.log(1.0-alpha);
        } else {
            
            // Contribution from start of sequence up to first recomb region
            if (recombs.get(1).getStartLocus()>0) {
                thisLogP += Math.log(pStartCF)
                        + (recombs.get(1).getStartLocus()-1)*Math.log(1-alpha);
            }  else {
                thisLogP += Math.log(1.0-pStartCF)
                        - Math.log(alpha);
            }
            
            // Contribution from remaining recomb regions and adjacent CF regions
            for (int ridx=1; ridx<recombs.size(); ridx++) {
                Recombination recomb = recombs.get(ridx);
                
                thisLogP += Math.log(alpha)
                        + (recomb.getEndLocus() - recomb.getStartLocus())*Math.log(1.0-deltaInv);
                
                if (ridx<recombs.size()-1) {
                    thisLogP += Math.log(deltaInv)
                            + (recombs.get(ridx+1).getStartLocus()-recomb.getEndLocus()-2)
                            *Math.log(1.0-alpha);
                } else {
                    if (recomb.getEndLocus()<sequenceLength-1) {
                        thisLogP += Math.log(deltaInv)
                                + (sequenceLength-1-recomb.getEndLocus()-1)
                                *Math.log(1.0-alpha);
                    }
                }
            }
        }
        
        
        return thisLogP;
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }
}
