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

import bacter.Recombination;
import beast.core.Description;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which reverses the nodes that an edge corresponds to "
        + "and the times at which it attaches to the clonal frame leaving "
        + "everything else unchanged.")
public class RecombinantEdgeFlip extends RecombinationGraphOperator {

    @Override
    public double proposal() {

        if (arg.getNRecombs()==0)
            return Double.NEGATIVE_INFINITY;
        
        Recombination recomb = arg.getRecombinations().get(
                Randomizer.nextInt(arg.getNRecombs())+1);
        
        Node node1 = recomb.getNode1();
        Node node2 = recomb.getNode2();
        
        double height1 = recomb.getHeight1();
        double height2 = recomb.getHeight2();

        if (node1 == node2 || node2.isRoot())
            return Double.NEGATIVE_INFINITY;
        
        if (height1<node2.getHeight()
                || height1>node2.getParent().getHeight()
                || height2<node1.getHeight()
                || height2>node1.getParent().getHeight())
            return Double.NEGATIVE_INFINITY;
        
        recomb.setNode1(node2);
        recomb.setNode2(node1);
        
        return 0.0;
    }
    
}
