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
import argbeast.RecombinationGraph;
import beast.core.Description;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Cause recombinant edge to hop between clonal frame edges.")
public class RecombinantEdgeHop extends RecombinationGraphOperator {

    @Override
    public double proposal() {

        RecombinationGraph arg = argInput.get();
        
        if (arg.getNRecombs()==0)
            return Double.NEGATIVE_INFINITY;
        
        // Select recombination at random
        Recombination recomb = arg.getRecombinations().get(
                Randomizer.nextInt(arg.getNRecombs())+1);
        
        // Select new attachment point:
        double u = Randomizer.nextDouble()*arg.getClonalFrameLength();
        Node nodeBelow;
        double newHeight;
        for (Node node : arg.getNodesAsArray()) {
            
        }
        
        return 0.0;
    }
    
}
