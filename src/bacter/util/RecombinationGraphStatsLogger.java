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

package bacter.util;

import bacter.Conversion;
import bacter.ConversionGraph;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import java.io.PrintStream;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class RecombinationGraphStatsLogger extends CalculationNode implements Loggable {

    public Input<ConversionGraph> argInput = new Input <ConversionGraph>(
            "arg", "Recombination graph to calculate summary statistics from.",
            Validate.REQUIRED);
    
    private ConversionGraph arg;
    
    @Override
    public void initAndValidate() {
        arg = argInput.get();
    }
    
    /**
     * Obtain mean length of converted regions described by ARG.
     * 
     * @param arg
     * @return mean length, or NaN if ARG has no recombinant edges.
     */
    public static double getMeanTractLength(ConversionGraph arg) {
        
        if (arg.getNConvs()<1)
            return Double.NaN;
        
        double mean = 0;
        for (Conversion recomb : arg.getConversions()) {
            if (recomb == null)
                continue;
            
            mean += recomb.getEndSite()-recomb.getStartSite()+1;
        }
        mean /= arg.getNConvs();
        
        return mean;
    }
    
    /**
     * Obtain mean number of loci between converted regions described
     * by ARG.
     * 
     * @param arg
     * @return mean count, or NaN if ARG has less than 2 conversions
     */
    public static double getMeanInterTractLength(ConversionGraph arg) {
        
        if (arg.getNConvs()<2)
            return Double.NaN;
        
        double mean = 0;
        for (int ridx=1; ridx<arg.getNConvs(); ridx++) {
            mean += arg.getConversions().get(ridx+1).getStartSite()
                    - arg.getConversions().get(ridx).getEndSite() - 1;
        }
        mean /= arg.getNConvs()-1;
        
        return mean;
    }
    
    /**
     * Obtain mean length of recombinant edges in ARG.
     * 
     * @param arg
     * @return mean length, or NaN if ARG has no recombinant edges
     */
    public static double getMeanEdgeLength(ConversionGraph arg) {
        
        if (arg.getNConvs()<1)
            return Double.NaN;
        
        double mean = 0.0;
        for (Conversion recomb : arg.getConversions()) {
            if (recomb == null)
                continue;
            
            mean += recomb.getHeight2()-recomb.getHeight1();
        }
        mean /= arg.getNConvs();
        
        return mean;
    }
    
    /**
     * Obtain mean height of point of departure of recombinant edges
     * in ARG.
     * 
     * @param arg
     * @return mean height, or NaN if ARG has no recombinant edges
     */
    public static double getMeanDepartureHeight(ConversionGraph arg) {
        if (arg.getNConvs()<1)
            return Double.NaN;
        
        double mean = 0.0;
        for (Conversion recomb : arg.getConversions()) {
            if (recomb == null)
                continue;
            
            mean += recomb.getHeight1();
        }
        mean /= arg.getNConvs();
        
        return mean;
    }
    
    @Override
    public void init(PrintStream out) throws Exception {
        String id = getID();
        if (id == null || id.matches("\\s*"))
            id = arg.getID();

        out.print(id + ".CFheight\t"
                + id + ".CFlength\t"
                + id + ".nRecomb\t"
                + id + ".meanTractLength\t"
                + id + ".meanInterTractLength\t"
                + id + ".meanEdgeLength\t"
                + id + ".meanDepartureHeight\t");
    }

    @Override
    public void log(int nSample, PrintStream out) {
        out.print(arg.getRoot().getHeight() + "\t"
                + arg.getClonalFrameLength() + "\t"
                + arg.getNConvs() + "\t"
                + RecombinationGraphStatsLogger.getMeanTractLength(arg) + "\t"
                + RecombinationGraphStatsLogger.getMeanInterTractLength(arg) + "\t"
                + RecombinationGraphStatsLogger.getMeanEdgeLength(arg) + "\t"
                + RecombinationGraphStatsLogger.getMeanDepartureHeight(arg) + "\t");
    }

    @Override
    public void close(PrintStream out) { }
    
}
