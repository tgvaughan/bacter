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
import beast.core.Description;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Cause recombinant edge to hop between clonal frame edges.")
public class ConvertedEdgeHopContemp extends ACGOperator {

    public ConvertedEdgeHopContemp() { }
    
    @Override
    public double proposal() {

        if (acg.getConvCount() == 0)
            return Double.NEGATIVE_INFINITY;

        // Select recombination at random
        Conversion conv = acg.getConversions().get(
                Randomizer.nextInt(acg.getConvCount()));

        // Choose whether to move departure or arrival point
        boolean moveDeparture = conv.getNode2().isRoot() || Randomizer.nextBoolean();

        double pointHeight = moveDeparture ? conv.getHeight1() : conv.getHeight2();
        Node convNode = moveDeparture ? conv.getNode1() : conv.getNode2();

        // Find list of CF edges alive at pointHeight
        List<Node> intersectingEdges = new ArrayList<>();
        for (Node node : acg.getNodesAsArray()) {
            if (node.isRoot()
                    || node == convNode
                    || node.getHeight() > pointHeight
                    || node.getParent().getHeight() < pointHeight) {
                continue;
            }

            intersectingEdges.add(node);
        }

        if (intersectingEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;

        // Select new attachment point:
        if (moveDeparture)
            conv.setNode1(intersectingEdges.get(Randomizer.nextInt(intersectingEdges.size())));
        else
            conv.setNode2(intersectingEdges.get(Randomizer.nextInt(intersectingEdges.size())));

        return 0.0;
    }
    
}
