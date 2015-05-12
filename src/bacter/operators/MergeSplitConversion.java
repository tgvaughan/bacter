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
package bacter.operators;

import bacter.Conversion;
import bacter.operators.ACGOperator;
import bacter.operators.EdgeCreationOperator;
import beast.evolution.alignment.Alignment;
import beast.util.Randomizer;

/**
 * Merge-split move for the unrestricted model.  The merge move selects
 * two conversions at random, deletes one of them and extends the converted
 * region of the other to include all of the sites (and more) of the other.
 * 
 * The operator only applies to conversions that attach to identical
 * pairs of clonal frame edges.
 *
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
public class MergeSplitConversion extends ACGOperator {

    @Override
    public double proposal() {

        Alignment alignment = chooseAlignment();

        if (Randomizer.nextBoolean())
            return mergeProposal(alignment);
        else
            return splitProposal(alignment);
    }

    /**
     * Perform merge portion of merge/split move.
     *
     * @param alignment alignment on which to apply move
     * @return log of move HR
     */
    private double mergeProposal(Alignment alignment) {

        double logHGF = 0.0;

        if (acg.getConvCount(alignment)<2)
            return Double.NEGATIVE_INFINITY;

        Conversion conv1 = acg.getConversions(alignment).get(
            Randomizer.nextInt(acg.getConvCount(alignment)));

        Conversion conv2;
        do {
            conv2 = acg.getConversions(alignment).get(
                Randomizer.nextInt(acg.getConvCount(alignment)));
        } while (conv2 == conv1);

        if (conv2.getNode1() != conv1.getNode1() || conv2.getNode2() != conv1.getNode2())
            return Double.NEGATIVE_INFINITY;

        logHGF -= Math.log(1.0/(acg.getConvCount(alignment)*(acg.getConvCount(alignment)-1)));

        int minStart = conv1.getStartSite() < conv2.getStartSite()
            ? conv1.getStartSite() : conv2.getStartSite();

        int maxEnd = conv1.getEndSite() > conv2.getEndSite()
            ? conv1.getEndSite() : conv2.getEndSite();

        logHGF += 2.0*Math.log(0.5/(maxEnd-minStart+1));

        logHGF += Math.log(1.0/conv1.getNode1().getLength());

        if (conv1.getNode2().isRoot()) {
            double lambda = 1.0/(conv1.getHeight2()-conv1.getNode2().getHeight());
            logHGF += -lambda*(conv2.getHeight2()-conv1.getNode2().getHeight())
                    + Math.log(lambda);
        } else {
            logHGF += Math.log(1.0/conv1.getNode2().getLength());
        }

        acg.deleteConversion(conv2);
        conv1.setStartSite(minStart);
        conv1.setEndSite(maxEnd);

        logHGF += Math.log(1.0/acg.getConvCount(alignment));

        return logHGF;
    }

    /**
     * Perform split component of merge/split move.
     *
     * @param alignment alignment on which to apply move
     * @return log of move HR
     */
    private double splitProposal(Alignment alignment) {

        double logHGF = 0.0;

        if (acg.getConvCount(alignment) == 0)
            return Double.NEGATIVE_INFINITY;

        Conversion conv1 = acg.getConversions(alignment).get(
            Randomizer.nextInt(acg.getConvCount(alignment)));

        logHGF -= Math.log(1.0/acg.getConvCount(alignment));

        Conversion conv2 = new Conversion();
        conv2.setNode1(conv1.getNode1());
        conv2.setNode2(conv1.getNode2());

        int s1, s2, e1, e2;

        int m1 = conv1.getStartSite() + Randomizer.nextInt(conv1.getSiteCount());
        int m2 = conv1.getStartSite() + Randomizer.nextInt(conv1.getSiteCount());
        
        if (Randomizer.nextBoolean()) {
            s1 = conv1.getStartSite();
            s2 = m1;
        } else {
            s1 = m1;
            s2 = conv1.getStartSite();
        }

        if (Randomizer.nextBoolean()) {
            e1 = conv1.getEndSite();
            e2 = m2;
        } else {
            e1 = m2;
            e2 = conv1.getEndSite();
        }

        if (e1<s1 || e2<s2)
            return Double.NEGATIVE_INFINITY;

        logHGF -= 2.0*Math.log(0.5/(conv1.getSiteCount()));

        conv1.setStartSite(s1);
        conv1.setEndSite(e1);

        conv2.setStartSite(s2);
        conv2.setEndSite(e2);

        conv2.setHeight1(conv1.getNode1().getHeight()
                + Randomizer.nextDouble()*conv1.getNode1().getLength());

        logHGF -= Math.log(1.0/conv1.getNode1().getLength());

        if (conv1.getNode2().isRoot()) {
            double lambda = 1.0/(conv1.getHeight2()-conv1.getNode2().getHeight());
            conv2.setHeight1(conv1.getNode2().getHeight()
                    + Randomizer.nextExponential(lambda));
            logHGF -= -lambda*(conv2.getHeight2()-conv2.getNode2().getHeight())
                    + Math.log(lambda);
        } else {
            conv2.setHeight2(conv1.getNode2().getHeight()
                    + Randomizer.nextDouble()*conv1.getNode2().getLength());
            logHGF -= Math.log(1.0/conv1.getNode2().getLength());
        }

        if (conv2.getHeight2() < conv2.getHeight1())
            return Double.NEGATIVE_INFINITY;

        acg.addConversion(conv2);

        logHGF += Math.log(1.0/(acg.getConvCount(alignment)*(acg.getConvCount(alignment)-1)));

        return logHGF;
    }
    
}
