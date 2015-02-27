/*
 * Copyright (C) 2015 Tim Vaughan (tgvaughan@gmail.com)
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
package bacter.operators.unrestricted;

import bacter.Conversion;
import bacter.operators.ACGOperator;
import beast.util.Randomizer;

/**
 * Merge-split move for the unrestricted model.  The merge move selects
 * two conversions at random, deletes one of them and extends the converted
 * region of the other to include all of the sites (and more) of the other.
 *
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
public class MergeSplitConversion extends ACGOperator {

    @Override
    public double proposal() {
        if (Randomizer.nextBoolean())
            return mergeProposal();
        else
            return splitProposal();
    }

    private double mergeProposal() {

        if (acg.getConvCount()<2)
            return Double.NEGATIVE_INFINITY;

        Conversion conv1 = acg.getConversions().get(
            Randomizer.nextInt(acg.getConvCount()));

        Conversion conv2;
        do {
            conv2 = acg.getConversions().get(
                Randomizer.nextInt(acg.getConvCount()));
        } while (conv2 == conv1);

        if (conv2.getNode1() != conv1.getNode1() || conv2.getNode2() != conv2.getNode2())
            return Double.NEGATIVE_INFINITY;

        int minStart = conv1.getStartSite() < conv2.getStartSite()
            ? conv1.getStartSite() : conv2.getStartSite();

        int maxEnd = conv1.getEndSite() > conv2.getEndSite()
            ? conv1.getEndSite() : conv2.getEndSite();

        acg.deleteConversion(conv2);
        conv1.setStartSite(minStart);
        conv1.setStartSite(maxEnd);

        return 0.0;
    }

    private double splitProposal() {

        if (acg.getConvCount()<1)
            return Double.NEGATIVE_INFINITY;

        Conversion conv1 = acg.getConversions().get(
            Randomizer.nextInt(acg.getConvCount()));

        Conversion conv2 = new Conversion();

        int s1, s2, e1, e2;

        

        return 0.0;
    }
    
}
