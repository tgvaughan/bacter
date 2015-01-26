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
import bacter.Region;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import java.io.PrintStream;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ConversionGraphStatsLogger extends CalculationNode implements Loggable {

    public Input<ConversionGraph> acgInput = new Input <>(
            "acg", "Conversion graph to calculate summary statistics from.",
            Validate.REQUIRED);
    
    private ConversionGraph acg;
    
    @Override
    public void initAndValidate() {
        acg = acgInput.get();
    }
    
    /**
     * Obtain mean length of converted regions described by ACG.
     * 
     * @param acg
     * @return mean length, or NaN if ACG has no converted edges.
     */
    public static double getMeanTractLength(ConversionGraph acg) {
        
        if (acg.getConvCount()<1)
            return Double.NaN;
        
        double mean = 0;
        for (Conversion conv : acg.getConversions())
            mean += conv.getEndSite()-conv.getStartSite()+1;

        mean /= acg.getConvCount();
        
        return mean;
    }

    /**
     * Obtain mean length of contiguous regions having the same marginal
     * tree.
     * 
     * @param acg
     * @return mean length
     */
    public static double getMeanRegionLength(ConversionGraph acg) {
        double sum = 0.0;
        for (Region region : acg.getRegions())
            sum += region.getRegionLength();

        return sum/acg.getRegionCount();
    }

    /**
     * Obtain average position of conversion starts.
     * 
     * @param acg
     * @return average site index or NaN if no conversions in ACG.
     */
    public static double getMeanStartSite(ConversionGraph acg) {
        if (acg.getConvCount()<1)
            return Double.NaN;

        double mean = 0.0;
        for (Conversion conv : acg.getConversions())
            mean += conv.getStartSite();

        mean /= acg.getConvCount();

        return mean;
    }

    /**
     * Obtain average position of conversion ends.
     * 
     * @param acg
     * @return average site index or NaN if no conversions in ACG.
     */
    public static double getMeanEndSite(ConversionGraph acg) {
        if (acg.getConvCount()<1)
            return Double.NaN;

        double mean = 0.0;
        for (Conversion conv : acg.getConversions())
            mean += conv.getEndSite();

        mean /= acg.getConvCount();

        return mean;
    }

    /**
     * Obtain mean length of converted edges in ACG.
     * 
     * @param acg
     * @return mean length, or NaN if ACG has no converted edges
     */
    public static double getMeanEdgeLength(ConversionGraph acg) {
        
        if (acg.getConvCount()<1)
            return Double.NaN;
        
        double mean = 0.0;
        for (Conversion conv : acg.getConversions())
            mean += conv.getHeight2()-conv.getHeight1();

        mean /= acg.getConvCount();
        
        return mean;
    }
    
    /**
     * Obtain mean height of point of departure of converted edges
     * in ACG.
     * 
     * @param acg
     * @return mean height, or NaN if ACG has no converted edges
     */
    public static double getMeanDepartureHeight(ConversionGraph acg) {
        if (acg.getConvCount()<1)
            return Double.NaN;
        
        double mean = 0.0;
        for (Conversion conv : acg.getConversions()) {
            mean += conv.getHeight1();
        }
        mean /= acg.getConvCount();
        
        return mean;
    }
    
    @Override
    public void init(PrintStream out) throws Exception {
        String id = getID();
        if (id == null || id.matches("\\s*"))
            id = acg.getID();

        out.print(id + ".CFheight\t"
                + id + ".CFlength\t"
                + id + ".nRecomb\t"
                + id + ".meanTractLength\t"
                + id + ".meanRegionLength\t"
                + id + ".meanStartSite\t"
                + id + ".meanEndSite\t"
                + id + ".meanEdgeLength\t"
                + id + ".meanDepartureHeight\t");
    }

    @Override
    public void log(int nSample, PrintStream out) {
        out.print(acg.getRoot().getHeight() + "\t"
                + acg.getClonalFrameLength() + "\t"
                + acg.getConvCount() + "\t"
                + ConversionGraphStatsLogger.getMeanTractLength(acg) + "\t"
                + ConversionGraphStatsLogger.getMeanRegionLength(acg) + "\t"
                + ConversionGraphStatsLogger.getMeanStartSite(acg) + "\t"
                + ConversionGraphStatsLogger.getMeanEndSite(acg) + "\t"
                + ConversionGraphStatsLogger.getMeanEdgeLength(acg) + "\t"
                + ConversionGraphStatsLogger.getMeanDepartureHeight(acg) + "\t");
    }

    @Override
    public void close(PrintStream out) { }
    
}
