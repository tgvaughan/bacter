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

package bacter;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class RegionList {

    private final List<Region> regions;
    private boolean dirty;

    /**
     * Ancestral conversion graph this list belongs to.
     */
    private final ConversionGraph acg;

    /**
     * Construct a new region list for the given ACG.  There should only
     * be one of these objects per ACG object, created during the ACG
     * initAndValidate().
     * 
     * @param acg Conversion graph from which to compute region list.
     */
    public RegionList(ConversionGraph acg) {
        this.acg = acg;
        regions = new ArrayList<>();
        dirty = true;
    }

    /**
     * Obtain list of contiguous regions having fixed marginal trees.
     * 
     * @return region list
     */
    public List<Region> getRegions() {
        updateRegionList();

        return regions;
    }

    /**
     * Mark the region list as dirty.
     */
    public void makeDirty() {
        dirty = true;
    }
   
    /**
     * Assemble list of regions of contiguous sites that possess a single
     * marginal tree.
     */
    public void updateRegionList() {
        if (!dirty)
            return;
        
        regions.clear();

        List<Conversion> convOrderedByStart = Lists.newArrayList();
        convOrderedByStart.addAll(acg.getConversions());
        convOrderedByStart.sort((Conversion o1, Conversion o2) -> o1.startSite - o2.startSite);

        List<Conversion> convOrderedByEnd = Lists.newArrayList();
        convOrderedByEnd.addAll(acg.getConversions());
        convOrderedByEnd.sort((Conversion o1, Conversion o2) -> o1.endSite - o2.endSite);

        Set <Conversion> activeConversions = Sets.newHashSet();

        Region currentRegion = new Region();
        currentRegion.leftBoundary = 0;

        while (!convOrderedByStart.isEmpty() || !convOrderedByEnd.isEmpty()) {

            int nextConvStartBoundary = Integer.MAX_VALUE;
            if (!convOrderedByStart.isEmpty())
                nextConvStartBoundary = convOrderedByStart.get(0).startSite;

            int nextConvEndBoundary = -1;
            if (!convOrderedByEnd.isEmpty())
                nextConvEndBoundary = convOrderedByEnd.get(0).endSite+1;

            if (nextConvStartBoundary<nextConvEndBoundary) {
                if (nextConvStartBoundary != currentRegion.leftBoundary) {
                    currentRegion.rightBoundary = nextConvStartBoundary;
                    regions.add(currentRegion);
                    currentRegion = new Region();
                    currentRegion.leftBoundary = nextConvStartBoundary;
                }

                activeConversions.add(convOrderedByStart.get(0));
                currentRegion.activeConversions.add(convOrderedByStart.get(0));
                convOrderedByStart.remove(0);

            } else {
                if (nextConvEndBoundary != currentRegion.leftBoundary) {
                    currentRegion.rightBoundary = nextConvEndBoundary;
                    regions.add(currentRegion);
                    currentRegion = new Region();
                    currentRegion.leftBoundary = nextConvEndBoundary;
                }

                activeConversions.remove(convOrderedByEnd.get(0));
                currentRegion.activeConversions.remove(convOrderedByEnd.get(0));
                convOrderedByEnd.remove(0);
            }
        }

        if (currentRegion.leftBoundary<acg.getSequenceLength()) {
            if (!currentRegion.isClonalFrame())
                throw new RuntimeException("Error updating region list!");

            currentRegion.rightBoundary = acg.getSequenceLength();
            regions.add(currentRegion);
        }

        dirty = false;
    } 
    
}
