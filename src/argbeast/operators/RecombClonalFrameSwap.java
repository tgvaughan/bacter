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
import com.google.common.collect.Lists;
import java.util.Collections;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which swaps the role of the clonal frame with that of "
        + "one of the marginal trees resulting from a conversion.  Note that "
        + "this move conserves the marginal trees themselves.")
public class RecombClonalFrameSwap extends EdgeCreationOperator {

    @Override
    public double proposal() {
        
        double logP = 0.0;

        if (arg.getNRecombs()==0)
            return Double.NEGATIVE_INFINITY;
        
        // 1. Choose recombination
        
        Recombination recomb = arg.getRecombinations().get(
                Randomizer.nextInt(arg.getNRecombs())+1);

        if (recomb.getNode1()==recomb.getNode2())
            return Double.NEGATIVE_INFINITY;
        
        logP -= Math.log(1.0/arg.getNRecombs());

        String oldARG = arg.getExtendedNewick(true);
        
        // 2. Make marginal tree of chosen recomb the new clonal frame.
        
        Node floating = recomb.getNode1().getParent();
        Node sister = floating.getLeft()==recomb.getNode1()
                ? floating.getRight()
                : floating.getLeft();
        
        floating.removeChild(sister);
        
        if (recomb.getNode2() == floating)
            recomb.setNode2(sister);
        
        if (floating.isRoot())
            sister.setParent(null);
        else {
            Node grandparent = floating.getParent();
            floating.setParent(null);
            grandparent.removeChild(floating);
            grandparent.addChild(sister);
        }
        
        Node newSister = recomb.getNode2();

        if (!newSister.isRoot()) {
            Node newParent = newSister.getParent();
            newParent.removeChild(newSister);
            newParent.addChild(floating);
        }
        floating.addChild(newSister);
        floating.setHeight(recomb.getHeight2());
        
        if (floating.isRoot())
            arg.setRoot(floating);
        
        // 3. Record the site ranges the existing conversions apply to and
        // delete those conversions. (The CF now applies to those conversions.)
        
        List<Integer> startSites = Lists.newArrayList();
        List<Integer> endSites = Lists.newArrayList();

        while (arg.getRecombinations().size()>1) {
            Recombination thisRecomb = arg.getRecombinations().get(arg.getNRecombs());
            startSites.add(thisRecomb.getStartSite());
            endSites.add(thisRecomb.getEndSite());
            arg.deleteRecombination(arg.getRecombinations().get(arg.getNRecombs()));
        }
        
        Collections.reverse(startSites);
        Collections.reverse(endSites);
        
        // 4. Create new conversions corresponding to the regions originally
        // governed by the CF.
        
        if (startSites.get(0)>0) {
            Recombination newRecomb = new Recombination();
            newRecomb.setStartSite(0);
            newRecomb.setEndSite(startSites.get(0)-1);
            logP -= attachEdge(newRecomb);
            arg.addRecombination(newRecomb);
        }
        
        for (int i=0; i<startSites.size()-1; i++) {
            Recombination newRecomb = new Recombination();
            newRecomb.setEndSite(endSites.get(i)+1);
            newRecomb.setEndSite(startSites.get(i+1)-1);
            logP -= attachEdge(newRecomb);
            arg.addRecombination(newRecomb);
        }
        
        if (endSites.get(endSites.size()-1)<arg.getSequenceLength()-1) {
            Recombination newRecomb = new Recombination();
            newRecomb.setStartSite(endSites.get(endSites.size()-1)+1);
            newRecomb.setEndSite(arg.getSequenceLength()-1);
            logP -= attachEdge(newRecomb);
            arg.addRecombination(newRecomb);
        }
        
        System.out.println(oldARG);
        System.out.println(arg.getExtendedNewick(true));
        
        return 0.0;
    }
    
}
