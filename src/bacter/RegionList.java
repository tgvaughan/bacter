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

import bacter.util.IntRanges;
import beast.evolution.tree.Node;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

/**
 * This class is used to maintain a list of marginal tree regions
 * corresponding to a given ACG.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class RegionList {

    private final List<Region> regions;
    private Locus locus;
    private boolean dirty;

    private Map<Conversion, List<Integer>> affectedSites = new HashMap<>();
    private Map<Conversion, Integer> affectedSiteCount = new HashMap<>();
    private Map<Conversion, Double> affectedSiteFraction = new HashMap<>();
    
    private ACGEventList acgEventList;
    
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
		synchronized (this) {
			if (dirty) {
				updateRegionList();
				dirty = false;
			}
		}

        return regions;
    }

    /**
     * Retrieve number of regions in region list.
     * 
     * @return region count
     */
//    public int getRegionCount() {
//        updateRegionList();
//
//        return regions.size();
//    }

    /**
     * Mark the region list as dirty.
     */
    public void makeDirty() {
        dirty = true;
    }
   
    
    public Map<Conversion, Integer> getAffectedSiteCount(){  
    	return null;
    }
    
    /**
     * Assemble list of regions of contiguous sites that possess a single
     * marginal tree.
     */
    private void updateRegionList() {
    	
		regions.clear();

		/*
		 * Assemble lists of conversions ordered by start and end sites. Note
		 * that these are COPIES of the conversion objects attached to the ACG.
		 * This ensures that subsequent modifications of these objects won't
		 * break our contract with the HashSet<Conversion> objects in the
		 * likelihood code.
		 */
		List<Conversion> convOrderedByStart = new ArrayList<>();
		Map<Conversion, Integer> affectedSiteCount = acg.affectedSiteList.getAffectedSiteCount();
		acg.getConversions(locus).forEach(conversion -> {
			if (affectedSiteCount.get(conversion) > 0)
				convOrderedByStart.add(conversion.getCopy());
		});
		convOrderedByStart.sort((Conversion o1, Conversion o2) -> o1.startSite - o2.startSite);

		List<Conversion> convOrderedByEnd = new ArrayList<>();
		convOrderedByEnd.addAll(convOrderedByStart);
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
				Region region = new Region(lastBoundary, nextBoundary, activeConversions);
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
			Region region = new Region(lastBoundary, locus.getSiteCount(), new HashSet<>());
			regions.add(region);
		}
	}
    
    private void updateAffectedSiteCounts(Locus locus){
    	
        affectedSites.clear();
        affectedSiteCount.clear();
        affectedSiteFraction.clear();
        
        Map<Node, List<Integer>> activeCFNodes = new HashMap<>();
        Set<Conversion> activeConversions = new HashSet<>();

        //List<Integer> ancestralSitesCF;

        int leavesSeen = 0;
        boolean mrcaReached = false;
        for (ACGEventList.Event event : acgEventList.getACGEvents()) {

            if (mrcaReached) {
                if (event.type == ACGEventList.EventType.CONV_DEPART) {
                    affectedSites.put(event.conversion, new ArrayList<>());
                    affectedSiteCount.put(event.conversion, 0);
                    affectedSiteFraction.put(event.conversion, 0.0);
                }

                continue;
            }

            switch (event.type) {
                case CF_LEAF:
                	activeCFNodes.put(event.node, getLeafAncestralSites());
                    leavesSeen += 1;
                    break;

                case CF_COALESCENCE:
                    Node node1 = event.node.getLeft();
                    Node node2 = event.node.getRight();

                    List<Integer> ancestralSitesCF = IntRanges.getUnion(activeCFNodes.get(node1), activeCFNodes.get(node2));                    

                    activeCFNodes.remove(node1);
                    activeCFNodes.remove(node2);
                    activeCFNodes.put(event.node, ancestralSitesCF);

                    if (leavesSeen == acg.getLeafNodeCount() && haveReachedAllMRCAs(activeCFNodes, activeConversions))
                        mrcaReached = true;

                    break;

                case CONV_DEPART:
                	if(event.conversion.getLocus() == locus){
                    List<Integer> inside = new ArrayList<>();
                    List<Integer> outside = new ArrayList<>();
                    IntRanges.partitionRanges(activeCFNodes.get(event.node),
                            event.conversion.getStartSite(),
                            event.conversion.getEndSite() + 1,
                            inside, outside);

                    affectedSites.put(event.conversion, inside);
                    affectedSiteCount.put(event.conversion, IntRanges.getTotalSites(inside));
                    affectedSiteFraction.put(event.conversion, IntRanges.getTotalSites(inside) / (double) event.conversion.getSiteCount());                
                    activeCFNodes.put(event.node, outside);
                    activeConversions.add(event.conversion);
                	}
                    break;

                case CONV_ARRIVE:
                	activeCFNodes.put(event.node, IntRanges.getUnion(affectedSites.get(event.conversion), activeCFNodes.get(event.node)));
                    activeConversions.remove(event.conversion);

                    if (leavesSeen == acg.getLeafNodeCount() && haveReachedAllMRCAs(activeCFNodes, activeConversions))
                        mrcaReached = true;
                    break;
            }

        }
    }
    
    /**
     * Assembles complete site list for association with a leaf node.
     *
     * @return list of sites
     */
    protected List<Integer> getLeafAncestralSites() {
    	List<Integer> siteRange = new ArrayList<>(2);
    	siteRange.add(0);
    	siteRange.add(locus.getSiteCount() - 1);
        return siteRange;
    }

    /**
     * Test to see whether MRCA of every site has been found.
     * This is actually pretty expensive.  There's got to be a better way...
     *
     * @param activeCFNodes set of active CF nodes and the sites they represent
     * @param activeConversions set of active conversions and the sites they represent
     * @return true if all sites have found an MRCA, false otherwise
     */
    protected boolean haveReachedAllMRCAs(Map<Node, List<Integer>> activeCFNodes, Set<Conversion> activeConversions) {

            List<Integer> startSites = new ArrayList<>();
            List<Integer> endSites = new ArrayList<>();
            for (Node node : activeCFNodes.keySet()) {
                for (int i = 0; i < activeCFNodes.get(node).size(); i += 2) {
                    startSites.add(activeCFNodes.get(node).get(i));
                    endSites.add(activeCFNodes.get(node).get(i + 1));
                }
            }

            for (Conversion conv : activeConversions) {
                for (int i = 0; i < affectedSites.get(conv).size(); i += 2) {
                    startSites.add(affectedSites.get(conv).get(i));
                    endSites.add(affectedSites.get(conv).get(i + 1));
                }
            }

            Collections.sort(startSites);
            Collections.sort(endSites);

            for (int i = 0; i < startSites.size() - 1; i++) {
                if (startSites.get(i + 1) < endSites.get(i))
                    return false;
            }
        return true;
    }
}
