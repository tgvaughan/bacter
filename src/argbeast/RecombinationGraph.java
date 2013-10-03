/*
 * Copyright (C) 2013 Tim Vaughan <tgvaughan@gmail.com>
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
package argbeast;

import beast.core.Description;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import com.google.common.collect.Lists;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Recombination graph based around the clonal frame.")
public class RecombinationGraph extends Tree {
    
    /**
     * List of recombinations on graph.
     */
    protected List<Recombination> recombs;
    
    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
        
        recombs = new ArrayList<Recombination>();
    }
    
    /**
     * Add recombination to graph, ensuring recombination list
     * remains sorted.
     * 
     * @param recomb 
     */
    public void addRecombination(Recombination recomb) {
        int i;
        for (i=0; i<recombs.size(); i++)
            if (recomb.startLocus>recombs.get(i).startLocus)
                break;
        
        recombs.add(i, recomb);
    }
    
    /**
     * Remove recombination from graph.
     * 
     * @param recomb 
     */
    public void deleteRecombination(Recombination recomb) {
        recombs.remove(recomb);
    }
    
    /**
     * Retrieve list of recombinations.
     * 
     * @return List of recombinations.
     */
    public List<Recombination> getRecombinations() {
        return recombs;
    }

    /**
     * Obtain number of recombination events.
     * 
     * @return Number of recombinations.
     */
    public int getNRecombs() {
        return recombs.size();
    }

    /**
     * Retrieve root for marginal tree defined by recomb.
     * 
     * @param recomb
     * @return root node
     */
    public Node getMarginalRoot(Recombination recomb) {
        
        if (recomb.node2.isRoot())
            return recomb.node1.getParent();
        else
            return getRoot();

    }
    
    /**
     * Retrieve parent of node in marginal tree defined by recomb.
     * 
     * @param node
     * @param recomb
     * @return node parent
     */
    public Node getMarginalParent(Node node, Recombination recomb) {
        
        if (node == recomb.node2) {
            return recomb.node1.getParent();
        }
        
        if (node.getParent() == recomb.node1.getParent())
            return node.getParent().getParent();
        
        return node.getParent();
    }
    
    /**
     * Retrieve children of node in marginal tree defined by recomb.
     * 
     * @param node
     * @param recomb
     * @return node children
     */
    public List<Node> getMarginalChildren(Node node, Recombination recomb) {
        List<Node> children = Lists.newArrayList();
        if (node == recomb.node1.getParent()) {
            children.add(recomb.node1);
            children.add(recomb.node2);
            return children;
        }
        
        for (Node child : node.getChildren()) {
            if (child == recomb.node2)
                children.add(recomb.node1.getParent());
            else
                children.add(child);
        }

        return children;
    }
}
