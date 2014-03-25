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

import beast.evolution.tree.Node;

/**
 * A class representing recombination events that are one-edge
 * modifications to the clonal frame.  A recombinant edge may leave and
 * connect with the same edge on the clonal frame, but no recombinant edges
 * may leave the root edge of the clonal frame.
 * 
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class Recombination {
    
    private RecombinationGraph arg;
    
    public Recombination() { }
    
    /**
     * Construct new recombination with specified properties.
     * @param node1
     * @param height1
     * @param node2
     * @param height2
     * @param startLocus
     * @param endLocus 
     */
    public Recombination(Node node1, double height1, Node node2, double height2,
            int startLocus, int endLocus) {
        this.node1 = node1;
        this.node2 = node2;
        this.height1 = height1;
        this.height2 = height2;
        this.startLocus = startLocus;
        this.endLocus = endLocus;
    }

    /**
     * Nodes below branches to which recombinant edge connects.
     */
    Node node1, node2;
    
    /**
     * Heights on branches at which recombinant edge connects.
     */
    double height1, height2;
    
    /**
     * Range of nucleotides affected by homologous gene conversion.
     */
    long startLocus, endLocus;
    
    /**
     * Obtain node below most recent point at which recombinant edge
     * attaches to clonal frame.
     * 
     * @return node
     */
    public Node getNode1() {
        return node1;
    }

    /**
     * Obtain node below oldest point at which recombinant edge attaches
     * to clonal frame.
     * 
     * @return node
     */
    public Node getNode2() {
        return node2;
    }
    
    /**
     * Return height at most recent attachment of recombinant edge to
     * clonal frame.
     * 
     * @return height
     */
    public double getHeight1() {
        return height1;
    }
    
    /**
     * Return height at oldest attachment of recombinant edge to clonal
     * frame.
     * @return height
     */
    public double getHeight2() {
        return height2;
    }
    
    /**
     * Return locus of start of alignment region affected by recombination.
     * @return locus
     */
    public long getStartLocus() {
        return startLocus;
    }
    
    /**
     * Return locus of end of alignment region affected by recombination.
     * @return locus
     */
    public long getEndLocus() {
        return endLocus;
    }

    /**
     * Set node below most recent point at which recombinant edge attaches
     * to clonal frame.
     * 
     * @param node1 
     */
    public void setNode1(Node node1) {
        startEditing();
        this.node1 = node1;
    }

    /**
     * Set node below oldest point at which recombinant edge attaches to
     * clonal frame.
     * 
     * @param node2 
     */
    public void setNode2(Node node2) {
        startEditing();
        this.node2 = node2;
    }

    /**
     * Set height at most recent attachment of recombinant edge to clonal
     * frame.
     * 
     * @param height1 
     */
    public void setHeight1(double height1) {
        startEditing();
        this.height1 = height1;
    }

    /**
     * Set height at oldest attachment of recombinant edge to clonal frame.
     * 
     * @param height2 
     */
    public void setHeight2(double height2) {
        startEditing();
        this.height2 = height2;
    }

    /**
     * Set locus of start of alignment region affected by conversion event.
     * 
     * @param startLocus 
     */
    public void setStartLocus(long startLocus) {
        startEditing();
        this.startLocus = startLocus;
    }

    /**
     * Set locus of end of alignment region affected by conversion event.
     * 
     * @param endLocus 
     */
    public void setEndLocus(long endLocus) {
        startEditing();
        this.endLocus = endLocus;
    }
    
    /**
     * Check validity of recombination specification: whether specified heights
     * belong to edges above specified nodes.
     * 
     * @return true if specification is valid
     */
    public boolean isValid() {
        if (height1>height2)
            return false;
        
        if (node1.getHeight()>height1)
            return false;
        
        if (node1.isRoot())
            return false;
        
        if (node1.getParent().getHeight()<height1)
            return false;
        
        if (node2.getHeight()>height2)
            return false;
        
        if (!node2.isRoot() && node2.getParent().getHeight()<height2)
            return false;
        
        if (startLocus>endLocus)
            return false;
        
        return true;
    }
    
    /**
     * Assign recombination graph.
     * @param arg 
     */
    public void setRecombinationGraph(RecombinationGraph arg) {
        this.arg = arg;
    }
    
    /**
     * Mark ARG statenode as dirty if available.
     */
    public void startEditing() {
        if (arg != null)
            arg.startEditing();
    }
    
    /**
     * Obtain new recombination with exactly the same
     * field values as this one.
     *
     * @return copy of Recombination object
     */
    public Recombination getCopy() {
        Recombination copy = new Recombination();
        copy.arg = arg;
        copy.startLocus = startLocus;
        copy.endLocus = endLocus;
        copy.node1 = node1;
        copy.node2 = node2;
        copy.height1 = height1;
        copy.height2 = height2;
        
        return copy;
    }

    /**
     * Returns true if obj is a Recombination object that represents an
     * identical recombination to this.
     * 
     * @param obj
     * @return true if recombinations are equivalent
     */
    @Override
    public boolean equals(Object obj) {
        if (obj instanceof Recombination) {
            Recombination otherRecomb = (Recombination)obj;
            return arg == otherRecomb.arg
                && startLocus == otherRecomb.startLocus
                && endLocus == otherRecomb.endLocus
                && node1 == otherRecomb.node1
                && node2 == otherRecomb.node2
                && height1 == otherRecomb.height1
                && height2 == otherRecomb.height2;
        } else
            return false;
    }

    /**
     * hashCode compatible with equals()
     *
     * @return hash
     */
    @Override
    public int hashCode() {
        int hash = 7;
        hash = 29 * hash + (this.arg != null ? this.arg.hashCode() : 0);
        hash = 29 * hash + (this.node1 != null ? this.node1.hashCode() : 0);
        hash = 29 * hash + (this.node2 != null ? this.node2.hashCode() : 0);
        hash = 29 * hash + (int) (Double.doubleToLongBits(this.height1) ^ (Double.doubleToLongBits(this.height1) >>> 32));
        hash = 29 * hash + (int) (Double.doubleToLongBits(this.height2) ^ (Double.doubleToLongBits(this.height2) >>> 32));
        hash = 29 * hash + (int) (this.startLocus ^ (this.startLocus >>> 32));
        hash = 29 * hash + (int) (this.endLocus ^ (this.endLocus >>> 32));
        return hash;
    }
}
