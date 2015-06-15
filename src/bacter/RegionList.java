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
    private Locus locus;
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
     * @param locus Locus with which this region list is associated.
     */
    public RegionList(ConversionGraph acg, Locus locus) {
        this.acg = acg;
        this.locus = locus;
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
     * Retrieve number of regions in region list.
     * 
     * @return region count
     */
    public int getRegionCount() {
        updateRegionList();

        return regions.size();
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

        List<Conversion> convOrderedByStart = new ArrayList<>();
        convOrderedByStart.addAll(acg.getConversions(locus));
        convOrderedByStart.sort((Conversion o1, Conversion o2) -> o1.startSite - o2.startSite);

        List<Conversion> convOrderedByEnd = new ArrayList<>();
        convOrderedByEnd.addAll(acg.getConversions(locus));
        convOrderedByEnd.sort((Conversion o1, Conversion o2) -> o1.endSite - o2.endSite);

        Set<Conversion> activeConversions = Sets.newHashSet();

        int lastBoundary = 0;

        while (!convOrderedByStart.isEmpty() || !convOrderedByEnd.isEmpty()) {

            int nextStart;
            if (!convOrderedByStart.isEmpty())
                nextStart = convOrderedByStart.get(0).getStartSite();
            else
                nextStart = Integer.MAX_VALUE;

            int nextEnd;
            if (!convOrderedByEnd.isEmpty())
                nextEnd = convOrderedByEnd.get(0).getEndSite() + 1;
            else
                nextEnd = Integer.MAX_VALUE;

            int nextBoundary = Math.min(nextStart, nextEnd);
            if (nextBoundary > lastBoundary) {
                Region region = new Region();
                region.leftBoundary = lastBoundary;
                region.rightBoundary = nextBoundary;
                region.activeConversions.addAll(activeConversions);
                regions.add(region);
            }

            if (nextStart < nextEnd) {
                activeConversions.add(convOrderedByStart.get(0));
                convOrderedByStart.remove(0);
                lastBoundary = nextStart;
            } else {
                activeConversions.remove(convOrderedByEnd.get(0));
                convOrderedByEnd.remove(0);
                lastBoundary = nextEnd;
            }
        }

        if (lastBoundary < locus.getSiteCount()) {
            Region region = new Region();
            region.leftBoundary = lastBoundary;
            region.rightBoundary = locus.getSiteCount();
            regions.add(region);
        }

        dirty = false;
    } 
    
}
