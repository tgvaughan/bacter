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

package bacter;

import beast.base.evolution.tree.Node;

/**
 * Adds a field to the Node class specifying whether a given
 * node in a marginal tree corresponds to a node in the clonal frame.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class MarginalNode extends Node {

    /**
     * If >= 0, denotes the node Nr of the CF node that
     * this node corresponds to.  If <0, this node does
     * not correspond to any CF node.
     */
    public int cfNodeNr = -1;
}
