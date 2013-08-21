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
import beast.evolution.tree.Tree;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Recombination graph based around the clonal frame.")
public class RecombinationGraph extends Tree {

    /**
     * Recombinations on graph.
     */
    protected List<Recombination> recombs;
    
    /**
     * List of marginal trees.
     */
    protected List<Tree> marginalTreeList;
    
    /**
     * Retrieve list of marginal trees.
     * @return 
     */
    public List<Tree> getMarginalTreeList() {

    }
    
}
