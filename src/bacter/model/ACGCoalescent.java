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
package bacter.model;

import bacter.CFEventList;
import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.Locus;
import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.coalescent.PopulationFunction;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.PoissonDistributionImpl;

import java.util.List;
import java.util.Random;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Appoximation to the coalescent with gene conversion.")
public class ACGCoalescent extends TreeDistribution {

    public Input<PopulationFunction> popFuncInput = new Input<>(
            "populationModel", "Population model.", Input.Validate.REQUIRED);
    
    public Input<RealParameter> rhoInput = new Input<>("rho",
            "Recombination rate parameter.", Input.Validate.REQUIRED);
    
    public Input<RealParameter> deltaInput = new Input<>("delta",
            "Tract length parameter.", Input.Validate.REQUIRED);

    public Input<Integer> lowerCCBoundInput = new Input<>("lowerConvCountBound",
            "Lower bound on conversion count.", 0);

    public Input<Integer> upperCCBoundInput = new Input<>("upperConvCountBound",
            "Upper bound on conversion count.", Integer.MAX_VALUE);

    // WARNING: The presence of the following input is a hack. This value is
    // not directly used by BEAST, but is instead copied via a custom BEAUti
    // connector to the similarly-named input to ConversionGraph.
    //
    // The reason for this ugliness is that BEAUti does not allow users to
    // alter inputs to Trees.  (The Tree input editor is used to set tip dates.)
    public Input<Boolean> wholeLocusConversionsInput = new Input<>(
            "wholeLocusConversionsOnly",
            "Only allow whole loci to be converted.", false);

    ConversionGraph acg;
    PopulationFunction popFunc;

    public ACGCoalescent() {
        treeInput.setRule(Input.Validate.REQUIRED);
    }
    
    @Override
    public void initAndValidate() {
        if (!(treeInput.get() instanceof ConversionGraph))
            throw new IllegalArgumentException("Tree input to ACGCoalescent " +
                    "must specify a ConversionGraph.");

        if (treeIntervalsInput.get() != null)
            throw new IllegalArgumentException("ACGCoalescent does not accept " +
                    "the treeIntervals input.");

        acg = (ConversionGraph)treeInput.get();
        popFunc = popFuncInput.get();
    }
    
    @Override
    public double calculateLogP() {

        // Check whether conversion count exceeds bounds.
        if (acg.getTotalConvCount()<lowerCCBoundInput.get()
                || acg.getTotalConvCount()>upperCCBoundInput.get())
            return Double.NEGATIVE_INFINITY;

        logP = calculateClonalFrameLogP();
        double poissonMean = rhoInput.get().getValue()
                *acg.getClonalFrameLength()
                *(acg.getTotalConvertibleSequenceLength()
                +acg.getConvertibleLoci().size()*(deltaInput.get().getValue()-1.0));

        // Probability of conversion count:
        if (poissonMean>0.0) {
            logP += -poissonMean + acg.getTotalConvCount()*Math.log(poissonMean);
            //      - GammaFunction.lnGamma(acg.getConvCount()+1);
        } else {
            if (acg.getTotalConvCount()>0)
                logP = Double.NEGATIVE_INFINITY;
        }
        

        for (Locus locus : acg.getConvertibleLoci())
            for (Conversion conv : acg.getConversions(locus))
                logP += calculateConversionLogP(conv);
        
        // This N! takes into account the permutation invariance of
        // the individual conversions, and cancels with the N! in the
        // denominator of the Poissonian above.
        // logP += GammaFunction.lnGamma(acg.getConvCount() + 1);

        if (lowerCCBoundInput.get()>0 || upperCCBoundInput.get()<Integer.MAX_VALUE) {
            try {
                logP -= new PoissonDistributionImpl(poissonMean)
                        .cumulativeProbability(
                                lowerCCBoundInput.get(),
                                upperCCBoundInput.get());
            } catch (MathException e) {
                throw new RuntimeException("Error computing modification to ARG " +
                        "prior density required by conversion number constraint.");
            }
        }

        return logP;
    }

    /**
     * Compute probability of clonal frame under coalescent.
     * 
     * @return log(P)
     */
    public double calculateClonalFrameLogP() {
        
        List<CFEventList.Event> events = acg.getCFEvents();
        
        double thisLogP = 0.0;
        
        for (int i=0; i<events.size()-1; i++) {
            double timeA = events.get(i).getHeight();
            double timeB = events.get(i+1).getHeight();

            double intervalArea = popFunc.getIntegral(timeA, timeB);
            int k = events.get(i).getLineageCount();
            thisLogP += -0.5*k*(k-1)*intervalArea;
            
            if (events.get(i+1).getType()==CFEventList.EventType.COALESCENCE)
                thisLogP += Math.log(1.0/popFunc.getPopSize(timeB));
        }
        
        return thisLogP;
    }
    
    /**
     * Compute probability of recombinant edges under conditional coalescent.
     * @param conv conversion with which edge is associated
     * @return log(P)
     */
    public double calculateConversionLogP(Conversion conv) {

        double thisLogP = 0.0;

        List<CFEventList.Event> events = acg.getCFEvents();

        // Probability density of location of recombinant edge start
        thisLogP += Math.log(1.0/acg.getClonalFrameLength());

        // Identify interval containing the start of the recombinant edge
        int startIdx = 0;
        while (events.get(startIdx+1).getHeight() < conv.getHeight1())
            startIdx += 1;

        for (int i=startIdx; i<events.size() && events.get(i).getHeight()<conv.getHeight2(); i++) {

            double timeA = Math.max(events.get(i).getHeight(), conv.getHeight1());

            double timeB;
            if (i<events.size()-1)
                timeB = Math.min(conv.getHeight2(), events.get(i+1).getHeight());
            else
                timeB = conv.getHeight2();

            double intervalArea = popFunc.getIntegral(timeA, timeB);
            thisLogP += -events.get(i).getLineageCount()*intervalArea;
        }

        // Probability of single coalescence event
        thisLogP += Math.log(1.0/popFunc.getPopSize(conv.getHeight2()));

        // Probability of start site:
        if (conv.getStartSite()==0) {
            thisLogP += Math.log(deltaInput.get().getValue()
                    / (acg.getConvertibleLoci().size() * (deltaInput.get().getValue() - 1)
                    + acg.getTotalConvertibleSequenceLength()));
        } else {
            if (!acg.wholeLocusModeOn())
                thisLogP += Math.log(
                        1.0 / (acg.getConvertibleLoci().size() * (deltaInput.get().getValue() - 1)
                                + acg.getTotalConvertibleSequenceLength()));
            else
                return Double.NEGATIVE_INFINITY;
        }

        // Probability of end site:
        if (conv.getEndSite() == conv.getLocus().getSiteCount()-1) {
            thisLogP += (conv.getLocus().getSiteCount()-1-conv.getStartSite())
                    *Math.log(1.0 - 1.0/deltaInput.get().getValue());
        } else {
            if  (!acg.wholeLocusModeOn())
                thisLogP += (conv.getEndSite()-conv.getStartSite())
                        *Math.log(1.0 - 1.0/deltaInput.get().getValue())
                        -Math.log(deltaInput.get().getValue());
            else
                return Double.NEGATIVE_INFINITY;
        }

        return thisLogP;
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    @Override
    public List<String> getArguments() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public List<String> getConditions() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
