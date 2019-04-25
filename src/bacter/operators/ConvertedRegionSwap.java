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

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Selects a pair of conversion events and swaps the recombinations"
        + "they correspond to.")
public class ConvertedRegionSwap extends ACGOperator {

    public ConvertedRegionSwap() { }
    
    @Override
    public double proposal() {
        
        if (acg.getTotalConvCount()<2)
            return Double.NEGATIVE_INFINITY;

        // Select a random pair of recombinations
        Conversion conv1 = chooseConversion();
        Conversion conv2;
        do {
            conv2 = chooseConversion();
        } while (conv2 == conv1);
        
        // Switch edges corresponding to recombinations.  (Can't switch
        // loci, as this would break ordering constraint.)
        
        double tmpHeight1 = conv1.getHeight1();
        double tmpHeight2 = conv1.getHeight2();
        Node tmpNode1 = conv1.getNode1();
        Node tmpNode2 = conv1.getNode2();
        
        conv1.setHeight1(conv2.getHeight1());
        conv1.setHeight2(conv2.getHeight2());
        conv1.setNode1(conv2.getNode1());
        conv1.setNode2(conv2.getNode2());
        
        conv2.setHeight1(tmpHeight1);
        conv2.setHeight2(tmpHeight2);
        conv2.setNode1(tmpNode1);
        conv2.setNode2(tmpNode2);

        assert !acg.isInvalid() : "CRS produced invalid state.";
        
        return 0.0;
    }
    
}
