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
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.alignment.Alignment;
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
     * Unlike Trees, Recombination graphs require an alignment to be specified
     * so that the regions of the alignment affected by recombination events
     * can be recorded.
     */
    public Input<Alignment> alignmentInput = new Input<Alignment>("alignment",
            "Sequence alignment corresponding to graph.", Validate.REQUIRED);
    
    protected int sequenceLength;
    
    /**
     * List of recombinations on graph.
     */
    protected List<Recombination> recombs;
    
    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
        
        recombs = new ArrayList<Recombination>();
        recombs.add(null); // Represents the clonal frame.
        
        sequenceLength = alignmentInput.get().getSiteCount();
    }
    
    /**
     * Retrieve length of sequence, identifying bounds of recombination loci.
     * 
     * @return sequence length
     */
    public int getSequenceLength() {
        return sequenceLength;
    }
    
    /**
     * Add recombination to graph, ensuring recombination list
     * remains sorted.
     * 
     * @param recomb 
     */
    public void addRecombination(Recombination recomb) {
        int i;
        for (i=1; i<recombs.size(); i++)
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
        if (recomb == null)
            throw new IllegalArgumentException("Cannot delete the clonal frame!");
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
        return recombs.size()-1;
    }
    
    /**
     * @return Total length of all edges in clonal frame.
     */
    public double getClonalFrameLength() {
        double length = 0.0;
        for (Node node : m_nodes) {
            if (node.isRoot())
                continue;
            length += node.getLength();
        }
        
        return length;
    }

    /**
     * Retrieve root for marginal tree defined by recomb.
     * 
     * @param recomb
     * @return root node
     */
    public Node getMarginalRoot(Recombination recomb) {
        if (recomb==null)
            return root;
        
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
        if (recomb==null)
            return node.getParent();
        
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
        if (recomb==null)
            return node.getChildren();
        
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

    /**
     * Height of node in marginal tree defined by recomb.
     * 
     * @param node
     * @param recomb
     * @return node height
     */
    public double getMarginalNodeHeight(Node node, Recombination recomb) {
        if (recomb==null)
            return node.getHeight();

        if (node == recomb.node1.getParent())
            return recomb.height2;
        else
            return node.getHeight();
    }
    
    /**
     * Length of edge between node and parent in marginal tree defined
     * by recomb.
     * 
     * @param node
     * @param recomb
     * @return edge length (zero if node is root of marginal tree)
     */
    public double getMarginalBranchLength(Node node, Recombination recomb) {
        if (recomb==null)
            return node.getLength();
        
        Node parent = getMarginalParent(node, recomb);
        if (parent == null)
            return 0.0;
        else
            return getMarginalNodeHeight(parent, recomb)
                    -getMarginalNodeHeight(node, recomb);
    }
    
    /**
     * Determine whether node is root of marginal tree defined by recomb.
     * 
     * @param node
     * @param recomb
     * @return true if node is marginal root, false otherwise.
     */
    public boolean isNodeMarginalRoot(Node node, Recombination recomb) {
        if (recomb==null)
            return node.isRoot();
        
        return (getMarginalParent(node, recomb) == null);
    }

    /**
     * Determine whether node is leaf in marginal tree defined by recomb.
     * Just here for aesthetic completeness - all marginal trees share the
     * same leaves.
     * 
     * @param node
     * @param recomb
     * @return true if node is leaf, false otherwise.
     */
    boolean isNodeMarginalLeaf(Node node, Recombination recomb) {
        return node.isLeaf();
    }
}
