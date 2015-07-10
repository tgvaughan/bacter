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
package bacter;

import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Node;
import java.util.Objects;

/**
 * A class representing recombination events that are one-edge
 * modifications to the clonal frame.  A recombinant edge may leave and
 * connect with the same edge on the clonal frame, but no recombinant edges
 * may leave the root edge of the clonal frame.
 * 
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class Conversion {

    /**
     * Conversion graph this node corresponds to.
     */
    public ConversionGraph acg;

    /**
     * Nodes below branches to which recombinant edge connects.
     */
    final public Node node1, node2;
    
    /**
     * Heights on branches at which recombinant edge connects.
     */
    final public double height1, height2;
    
    /**
     * Range of nucleotides affected by homologous gene conversion.
     */
    final public int startSite, endSite;

    /**
     * Locus with which conversion is associated.
     */
    final public Locus locus;

    /**
     * Used by ACGAnnotator to incoroporate additional metadata into
     * the summary ACG.
     */
    public String newickMetaDataBottom, newickMetaDataMiddle, newickMetaDataTop;

    /**
     * Construct new recombination with specified properties.
     *
     * @param node1
     * @param height1
     * @param node2
     * @param height2
     * @param startSite
     * @param endSite
     * @param locus
     */
    public Conversion(Node node1, double height1, Node node2, double height2,
            int startSite, int endSite, ConversionGraph acg, Locus locus) {
        this.node1 = node1;
        this.node2 = node2;
        this.height1 = height1;
        this.height2 = height2;
        this.startSite = startSite;
        this.endSite = endSite;
        this.locus = locus;
        this.acg = acg;
    }

    /**
     * @return total number of sites affected by this conversion.
     */
    public int getSiteCount() {
        return (int)(endSite - startSite + 1);
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
        
        if (startSite>endSite)
            return false;
        
        return true;
    }

    /**
     * Obtain new recombination with exactly the same
     * field values as this one.
     *
     * @return copy of Conversion object
     */
    public Conversion getCopy() {
        return new Conversion(
                node1,
                height1,
                node2,
                height2,
                startSite,
                endSite,
                acg,
                locus
        );
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Conversion that = (Conversion) o;

        if (Double.compare(that.height1, height1) != 0) return false;
        if (Double.compare(that.height2, height2) != 0) return false;
        if (startSite != that.startSite) return false;
        if (endSite != that.endSite) return false;
        if (!acg.equals(that.acg)) return false;
        if (!node1.equals(that.node1)) return false;
        if (!node2.equals(that.node2)) return false;
        return locus.equals(that.locus);

    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = acg.hashCode();
        result = 31 * result + node1.hashCode();
        result = 31 * result + node2.hashCode();
        temp = Double.doubleToLongBits(height1);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(height2);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        result = 31 * result + startSite;
        result = 31 * result + endSite;
        result = 31 * result + locus.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return String.format("Depart: (Node %d, height %g, site %d) "
                + "Arrive: (Node %d, height %g, site %d)",
                node1.getNr(), height1, startSite,
                node2.getNr(), height2, endSite);
    }
}
