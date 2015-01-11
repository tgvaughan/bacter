/*
 * Copyright (C) 2014 Tim Vaughan <tgvaughan@gmail.com>
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
import beast.core.Description;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Selects a pair of conversion events and swaps the recombinations"
        + "they correspond to.")
public class ConvertedRegionSwap extends RecombinationGraphOperator {

    public ConvertedRegionSwap() { }
    
    @Override
    public double proposal() {
        
        if (acg.getNConvs()<2)
            return Double.NEGATIVE_INFINITY;

        // Select a random pair of recombinations
        int[] idx = Randomizer.shuffled(acg.getNConvs());
        Conversion recomb1 = acg.getConversions().get(idx[0]+1);
        Conversion recomb2 = acg.getConversions().get(idx[1]+1);
        
        // Switch edges corresponding to recombinations.  (Can't switch
        // loci, as this would break ordering constraint.)
        
        double tmpHeight1 = recomb1.getHeight1();
        double tmpHeight2 = recomb1.getHeight2();
        Node tmpNode1 = recomb1.getNode1();
        Node tmpNode2 = recomb1.getNode2();
        
        recomb1.setHeight1(recomb2.getHeight1());
        recomb1.setHeight2(recomb2.getHeight2());
        recomb1.setNode1(recomb2.getNode1());
        recomb1.setNode2(recomb2.getNode2());
        
        recomb2.setHeight1(tmpHeight1);
        recomb2.setHeight2(tmpHeight2);
        recomb2.setNode1(tmpNode1);
        recomb2.setNode2(tmpNode2);
        
        return 0.0;
    }
    
}
