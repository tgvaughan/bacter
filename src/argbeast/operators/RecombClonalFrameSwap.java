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
        Node origNode2 = recomb.getNode2();
        double origHeight2 = recomb.getHeight2();
        double newHeight2 = origNode1Par.getHeight();

        String oldARG = arg.getExtendedNewick(true);
        
        for (Recombination otherRecomb : arg.getRecombinations()) {
            if (otherRecomb == null || otherRecomb == recomb)
                continue;

            // Can't perform move reversibly if there are recombinant edge
            // connections above height2 on node1
            if ((otherRecomb.getNode1() == origNode1
                    && otherRecomb.getHeight1()>recomb.getHeight2())
                    || (otherRecomb.getNode2() == origNode1
                    && otherRecomb.getHeight2()>recomb.getHeight2()))
                return Double.NEGATIVE_INFINITY;
            
            // Move all recombinant edge connections on origNode1Par to origNode1Sis
            
            if (otherRecomb.getNode1() == origNode1Par)
                otherRecomb.setNode1(origNode1Sis);
            
            if (otherRecomb.getNode2() == origNode1Par)
                otherRecomb.setNode2(origNode1Sis);
            
                    
            // Move all recombinant edge connetions on origNode2 above height2
            // to origNode1Par
            
            if (otherRecomb.getNode1() == origNode2 && otherRecomb.getHeight1()>recomb.getHeight2())
                otherRecomb.setNode1(origNode1Par);
            
            if (otherRecomb.getNode2() == origNode2 && otherRecomb.getHeight2()>recomb.getHeight2())
                otherRecomb.setNode2(origNode1Par);
        }
        
        // Ensure that if recomb attaches above node1Par it is reattached to node1Sis
        if (recomb.getNode2()==origNode1Par) {
            recomb.setNode2(origNode1Sis);
            origNode2 = origNode1Sis;
        }
        
        // Clonal frame topology changes
        
        origNode1Par.removeChild(origNode1Sis);
        origNode1Sis.setParent(origNode1Par.getParent());
        if (origNode1Sis.getParent() != null) {
            origNode1Sis.getParent().removeChild(origNode1Par);
            origNode1Sis.getParent().addChild(origNode1Sis);
        }

        origNode1Par.setParent(origNode2.getParent());
        origNode1Par.addChild(origNode2);
        if (origNode1Par.getParent() != null) {
            origNode1Par.getParent().removeChild(origNode2);
            origNode1Par.getParent().addChild(origNode1Par);
        }
        
        origNode1Par.setHeight(origHeight2);
        
        // Move original recombinant edge
        recomb.setHeight2(newHeight2);
        recomb.setNode2(origNode1Sis);
        
        System.out.println(oldARG);
        System.out.println(arg.getExtendedNewick(true));

        return 0.0;
    }
}
