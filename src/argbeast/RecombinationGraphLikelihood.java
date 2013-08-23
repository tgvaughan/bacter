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
import beast.evolution.likelihood.TreeLikelihood;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Probability of sequence data given recombination graph.")
public class RecombinationGraphLikelihood extends TreeLikelihood {

    List<int[]> patternWeights = new ArrayList<int[]>();
    int[] regionIndex;
    
    RecombinationGraph arg;
    
    @Override
    public void initAndValidate() throws Exception {
        if (!(m_tree.get() instanceof RecombinationGraph))
            throw new IllegalArgumentException("RecombinationGraphLikelihood "
                    + "can only be applied to RecombinationGraphs.");
        
        arg = (RecombinationGraph)m_tree.get();
        regionIndex = new int[m_data.get().getSiteCount()];
        
        calculatePatternWeights();
        
        super.initAndValidate();
    }
    
    @Override
    public double calculateLogP() {
        
        // Loop over sites
        
        return 0.0;
    }
    
    /**
     * Map each sight onto a region index. Sites marked with -1 belong to
     * the clonal frame.
     */
    private void calculateRegionMap() {
        
        int j=0;
        for (int ridx=0; ridx<arg.getRecombinations().size(); ridx++) {
            Recombination recomb = arg.getRecombinations().get(ridx);
            
            while (j < recomb.startLocus) {
                regionIndex[j] = -1;
                j += 1;
            }
            
            while (j <= recomb.endLocus) {
                regionIndex[j] = ridx;
            }
        }
    }

    /**
     * Calculate number of times each pattern coincides with each of
     * the potentially distinct genealogies.
     */
    private void calculatePatternWeights() {
        
        patternWeights.clear();
        for (Recombination recomb : arg.getRecombinations()) {
            m_data.get().getSiteCount();
        }
    }
}
