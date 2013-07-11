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

import beast.core.Plugin;
import beast.evolution.tree.Node;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class Recombination extends Plugin {

    /**
     * Nodes below branches to which recombinant edge connects.
     */
    protected Node node1, node2;
    
    /**
     * Heights on branches at which recombinant edge connects.
     */
    private double height1, height2;
    
    /**
     * Range of nucleotides affected by homologous gene conversion.
     */
    private int startLocus, endLocus;
    
}
