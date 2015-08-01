/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
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
import bacter.Locus;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

/**
 * Implementation of Wilson-Balding operator modified for the clonal frame
 * of the ACG. This version only proposes CF topology changes that are already
 * present in the form of conversions.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("CF/conversion swap operator.")
public class CFConversionSwap extends CFOperator {

    public Input<Boolean> includeRootInput = new Input<>("includeRoot",
            "Whether to include root variants of move.", true);

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
    }

    @Override
    public double proposal() {
        double logHGF = 0.0;

        // Determine whether we can apply this operator:
        if (acg.getLeafNodeCount()<3 || acg.getTotalConvCount()==0)
            return Double.NEGATIVE_INFINITY;

        // Acquire list of conversions compatible with swap:
        List<Conversion> compatible = getCompatibleConversions();
        if (compatible.isEmpty())
            return Double.NEGATIVE_INFINITY;

        // Select conversion to replace
        Conversion replaceConversion = compatible.get(
                Randomizer.nextInt(compatible.size()));
        logHGF -= Math.log(1.0/compatible.size());

        acg.deleteConversion(replaceConversion);

        Node srcNode = replaceConversion.getNode1();
        Node destNode = replaceConversion.getNode2();

        Node srcNodeP = srcNode.getParent();
        Node srcNodeS = getSibling(srcNode);
        double t_srcNodeP = srcNodeP.getHeight();

        if (destNode == srcNode.getParent())
            destNode = srcNodeS;

        double newTime = replaceConversion.getHeight2();

        // Conversion modification:
        replaceConversion.setNode2(srcNodeS);
        replaceConversion.setHeight2(t_srcNodeP);

        logHGF += getAffectedRegionProb(replaceConversion);
        logHGF -= drawAffectedRegion(replaceConversion);

        acg.addConversion(replaceConversion);

        if (!includeRootInput.get() && (srcNodeP.isRoot() || destNode.isRoot()))
            return Double.NEGATIVE_INFINITY;

        // Perform necessary conversion expansions/collapses:
        if (newTime < t_srcNodeP) {
            logHGF += collapseConversions(srcNode, destNode, newTime);
        } else {
            logHGF -= expandConversions(srcNode, destNode, newTime);
        }

        logHGF += Math.log(1.0/getCompatibleConversions().size());

        assert !acg.isInvalid() : "CFCS proposed invalid state.";

        return logHGF;
    }


    /**
     * @return list of conversions compatible with operator.
     */
    private List<Conversion> getCompatibleConversions() {
        List<Conversion> compatible = new ArrayList<>();
        for (Locus locus : acg.getLoci()) {
            for (Conversion conv : acg.getConversions(locus)) {
                if (conv.getNode1() != conv.getNode2())
                    compatible.add(conv);
            }
        }

        return compatible;
    }
}
