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
import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.math.GammaFunction;
import feast.input.In;
import java.util.List;
import java.util.Random;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Appoximation to the coalescent with gene conversion.")
public class ACGCoalescent extends ACGDistribution {
    
    public Input<ConversionGraph> acgInput = new In<ConversionGraph>(
            "acg", "Conversion graph.").setRequired();
    
    public Input<PopulationFunction> popFuncInput = new In<PopulationFunction>(
            "populationModel", "Population model.").setRequired();
    
    public Input<RealParameter> rhoInput = new In<RealParameter>("rho",
            "Recombination rate parameter.").setRequired();
    
    public Input<RealParameter> deltaInput = new In<RealParameter>("delta",
            "Tract length parameter.").setRequired();
    
    ConversionGraph acg;
    PopulationFunction popFunc;
    int sequenceLength;
    
    public ACGCoalescent() { }
    
    @Override
    public void initAndValidate() throws Exception {

        acg = acgInput.get();
        
        sequenceLength = acg.getSequenceLength();

        popFunc = popFuncInput.get();
        
        super.initAndValidate();
    }
    
    @Override
    public double calculateLogP() throws Exception {

        if (acg.isValid()) {

            logP = calculateClonalFrameLogP();

            // Probability of conversion count:
            double poissonMean = rhoInput.get().getValue()
                *acg.getClonalFrameLength()
                *(acg.getSequenceLength()+deltaInput.get().getValue());
            logP += -poissonMean + acg.getConvCount()*Math.log(poissonMean)
                - GammaFunction.lnGamma(acg.getConvCount()+1);

            for (Conversion conv : acg.getConversions())
                logP += calculateConversionLogP(conv);
            
        } else {
            logP = Double.NEGATIVE_INFINITY;
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
     * @param conv
     * @return log(P)
     */
    public double calculateConversionLogP(Conversion conv) {
        
        List<CFEventList.Event> events = acg.getCFEvents();
        
        // Probability density of location of recombinant edge start
        double thisLogP = Math.log(1.0/acg.getClonalFrameLength());

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
        if (conv.getStartSite()==0)
            thisLogP += Math.log((deltaInput.get().getValue() + 1)
                /(deltaInput.get().getValue() + acg.getSequenceLength()));
        else
            thisLogP += Math.log(1.0/(deltaInput.get().getValue()
                + acg.getSequenceLength()));

        // Probability of end site:
        int length = conv.getStartSite() - conv.getEndSite() + 1;
        double probEnd = Math.pow(1.0-1.0/deltaInput.get().getValue(), length-1)
            / deltaInput.get().getValue();
        
        // Include probability of going past the end:
        if (conv.getEndSite() == acg.getSequenceLength()-1)
            probEnd += Math.pow(1.0-1.0/deltaInput.get().getValue(),
                    acg.getSequenceLength()-conv.getStartSite());

        thisLogP += Math.log(probEnd);
        
        return thisLogP;
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    @Override
    public List<String> getArguments() {
        return null; // Doesn't seem to be used...
    }

    @Override
    public List<String> getConditions() {
        return null; // Doesn't seem to be used...
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
