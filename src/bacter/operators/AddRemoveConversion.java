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
package bacter.operators;

import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.util.RandomizedAlignment;
import beast.core.Description;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.util.Randomizer;

import java.io.PrintStream;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which adds and removes conversions to/from an ACG.")
public class AddRemoveConversion extends ConversionCreationOperator {
    
    public AddRemoveConversion() { }
    
    @Override
    public double proposal() {
        double logHGF = 0;

        Alignment alignment = chooseAlignment();
        
        if (Randomizer.nextBoolean()) {
            
            // Add
            
            logHGF += Math.log(1.0/(acg.getConvCount(alignment)+1));
            logHGF -= drawNewConversion(alignment);
            
        } else {
            
            // Remove
            
            if (acg.getConvCount(alignment)==0)
                return Double.NEGATIVE_INFINITY;
            
            // Select conversion to remove:
            Conversion conv = acg.getConversions(alignment).get(
                    Randomizer.nextInt(acg.getConvCount(alignment)));
            
            // Calculate HGF
            logHGF += getConversionProb(conv);
            logHGF -= Math.log(1.0/acg.getConvCount(alignment));
            
            // Remove conversion
            acg.deleteConversion(conv);

        }

        return logHGF;
    }
    
    /**
     * Add new conversion to ACG, returning the probability density of the
     * new edge and converted region location.
     *
     * @param alignment alignment with which to associate conversion
     * @return log of proposal density 
     */
    public double drawNewConversion(Alignment alignment) {
        Conversion newConversion = new Conversion(alignment);

        double logP = attachEdge(newConversion) + drawAffectedRegion(newConversion);

        acg.addConversion(newConversion);
        
        return logP;
    }
      
    /**
     * Obtain proposal density for the move which results in the current state
     * by adding the conversion conv to a state without that recombination.
     * 
     * @param conv conversion
     * @return log of proposal density
     */
    public double getConversionProb(Conversion conv) {
        return getEdgeAttachmentProb(conv) + getAffectedRegionProb(conv);
    }

    public static void main(String[] args) throws Exception {

        ConversionGraph acg = new ConversionGraph();
        ConstantPopulation popFunc = new ConstantPopulation();


        AddRemoveConversion operator = new AddRemoveConversion();
        operator.initByName("weight", 1.0,
            "acg", acg,
            "populationModel", popFunc,
            "rho", new RealParameter(Double.toString(1.0/10000.0)),
            "delta", new RealParameter("50.0"));
        popFunc.initByName("popSize", new RealParameter("1.0"));

        TaxonSet taxonSet = new TaxonSet();
        taxonSet.taxonsetInput.setValue(new Taxon("t1"), taxonSet);
        taxonSet.taxonsetInput.setValue(new Taxon("t2"), taxonSet);

        RandomizedAlignment alignment = new RandomizedAlignment();
        alignment.initByName("taxonSet", taxonSet, "sequenceLength", 10000);

        try (PrintStream ps = new PrintStream("out.txt")) {
            for (int i=0; i<100000; i++) {
                acg.initByName(
                    "alignment", alignment,
                    "fromString", "(0:1.0,1:1.0)2:0.0;");
                operator.drawNewConversion(alignment);
                
                ps.println(acg.getConversions(alignment).get(0).getStartSite() + " "
                    + acg.getConversions(alignment).get(0).getEndSite());
            }
        }
    }
}
