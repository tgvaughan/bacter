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

package argbeast.operators;

import argbeast.Recombination;
import beast.core.Description;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which swaps the role of the clonal frame with that of "
        + "one of the marginal trees resulting from a conversion.  Note that "
        + "this move conserves the marginal trees themselves.")
public class RecombClonalFrameSwap extends RecombinationGraphOperator {

    @Override
    public double proposal() {

        if (arg.getNRecombs()==0)
            return Double.NEGATIVE_INFINITY;
        
        Recombination recomb = arg.getRecombinations().get(
                Randomizer.nextInt(arg.getNRecombs())+1);

        if (recomb.getNode1()==recomb.getNode2())
            return Double.NEGATIVE_INFINITY;
        
        Node origNode1 = recomb.getNode1();
        Node origNode1Par = origNode1.getParent();
        Node origNode1Sis = getSibling(origNode1);
        
        return 0.0;
    }
}
