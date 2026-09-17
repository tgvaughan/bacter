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
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.TreeDistribution;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.util.GammaFunction;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.PoissonDistributionImpl;

import java.util.List;
import java.util.Random;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Approximation to the coalescent with gene conversion.")
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
    //circular genome mode edit
    public Input<Boolean> circularGenomeInput = new Input<>(
            "circularGenome",
            "The alignment is a circular genome", false);
    public Input<Boolean> betaBinomialEndSiteInput = new Input<>(
            "endSiteBetaBinom",
            "The prior for the end site of a conversion is a beta-binomial distribution.", false);

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

        // The following condition makes sure that in the case of a complete genome
        // the mean conversion length is smaller than half of the genome length
        // (following the convention of defining the shorter sequence part as conversion).
        //circular genome mode edit
        if (circularGenomeInput.get()){
            if (deltaInput.get().getValue() >= 0.5 * acg.getTotalConvertibleSequenceLength())
                throw new IllegalArgumentException("Delta prior input " +
                        "must be smaller than half of the genome length.");
            if (deltaInput.get().getUpper()>= 0.5 * acg.getTotalConvertibleSequenceLength()) {
                deltaInput.get().setUpper(Math.floor((acg.getTotalConvertibleSequenceLength() - 1.) * 0.5));
                System.out.println("Upper bound of delta is set to " + (Math.floor((acg.getTotalConvertibleSequenceLength() - 1.) * 0.5)));
            }
            if (!acg.circularGenomeModeOn()) {
                throw new IllegalArgumentException("Error: Circular genome mode turned on in ACGCoalescent but not in ConversionGraph (acg). Aborting. ");
            }
        }

        if (!circularGenomeInput.get() && acg.circularGenomeModeOn()) {
            throw new IllegalArgumentException("Error: Circular genome mode turned on in ConversionGraph (acg) but not in ACGCoalescent. Aborting. ");
        }

        //acg = (ConversionGraph)treeInput.get();
        popFunc = popFuncInput.get();
    }
    
    @Override
    public double calculateLogP() {

        // Check whether conversion count exceeds bounds.
        if (acg.getTotalConvCount()<lowerCCBoundInput.get()
                || acg.getTotalConvCount()>upperCCBoundInput.get())
            return Double.NEGATIVE_INFINITY;

        //circular genome mode edit
        logP = calculateClonalFrameLogP();
        double poissonMean = rhoInput.get().getValue()
                *acg.getClonalFrameLength()
                *(acg.getTotalConvertibleSequenceLength()
                + ( acg.circularGenomeModeOn() ? 0 :  acg.getConvertibleLoci().size()*(deltaInput.get().getValue()-1.0) )
        );

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
        //circular genome mode edit
        if (acg.circularGenomeModeOn()) {
            thisLogP += Math.log(1.0 / acg.getTotalConvertibleSequenceLength());
        } else if (conv.getStartSite()==0) {
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
        //circular genome mode edit
        if (acg.circularGenomeModeOn()) {
            if (acg.endSiteBetaBinomOn()) {
                int halfGenomeLength = (int) Math.floor((acg.getTotalConvertibleSequenceLength() - 1.) * 0.5);
                int kBetaBinom = conv.getSiteCount() - 1;
                double aBetaBinom = halfGenomeLength / (halfGenomeLength - deltaInput.get().getValue());
                double bBetaBinom = halfGenomeLength / deltaInput.get().getValue();
                thisLogP += GammaFunction.lnGamma(halfGenomeLength + 1) - GammaFunction.lnGamma(kBetaBinom + 1)
                        - GammaFunction.lnGamma(halfGenomeLength - kBetaBinom + 1)
                        + GammaFunction.lnGamma(kBetaBinom + aBetaBinom)
                        + GammaFunction.lnGamma(halfGenomeLength - kBetaBinom + bBetaBinom) - GammaFunction.lnGamma(bBetaBinom)
                        - GammaFunction.lnGamma(halfGenomeLength + aBetaBinom + bBetaBinom)
                        + GammaFunction.lnGamma(aBetaBinom + bBetaBinom) - GammaFunction.lnGamma(aBetaBinom);
            } else {
                thisLogP += (conv.getSiteCount() - 1)
                        *Math.log(1.0 - 1.0/deltaInput.get().getValue())
                        -Math.log(deltaInput.get().getValue())
                        -Math.log(1.0-Math.pow(1.0-1.0/deltaInput.get().getValue(), (int) Math.floor((acg.getTotalConvertibleSequenceLength()) * 0.5)));
            }
        } else if (conv.getEndSite() == conv.getLocus().getSiteCount()-1) {
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

    public static void main(String[] args) {

        int n = 3754;
        int k = 456;
        double delt = 410.0;
        double alpha = n/(n-delt);
        double beta = n/delt;

        long startTime = System.nanoTime();
        for(int j=0; j<1000000; j++) {
            double probConvGamma = GammaFunction.lnGamma(n + 1) - GammaFunction.lnGamma(k + 1) - GammaFunction.lnGamma(n - k + 1)
                    + GammaFunction.lnGamma(k + alpha) + GammaFunction.lnGamma(n - k + beta) - GammaFunction.lnGamma(beta)
                    - GammaFunction.lnGamma(n + alpha + beta) + GammaFunction.lnGamma(alpha + beta) - GammaFunction.lnGamma(alpha);
        }
        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
        System.out.println(duration/Math.pow(10,9));
    }
}
