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

import beast.base.evolution.tree.Node;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Maintains an ordered list of events which make up the clonal frame.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class CFEventList {

    /**
     * Class of events types on clonal frame.
     */
    public enum EventType {COALESCENCE, SAMPLE }

    /**
     * Class of events on clonal frame.
     */
    public static class Event {
        EventType type;
        double t;
        int lineages;
        Node node;

        /**
         * Construct event object corresponding to chosen CF node.
         *
         * @param node chosen CF node
         */
        public Event(Node node) {
            this.node = node;

            if (node.isLeaf())
                type = EventType.SAMPLE;
            else
                type = EventType.COALESCENCE;
            
            t = node.getHeight();
        }

        /**
         * Construct an event object with only the time defined.
         * Useful for conducting binary searches over the event list.
         *
         * @param t age of the event object
         */
        public Event(double t) {
            this.t = t;
        }

        public double getHeight() {
            return t;
        }
        
        public EventType getType() {
            return type;
        }

        public Node getNode() {
            return node;
        }
        
        /**
         * @return number of lineages _above_ this event.
         */
        public int getLineageCount() {
            return lineages;
        }
        
        @Override
        public String toString() {
            return "t: " + t + ", k: " + lineages + ", type: " + type;
        }
    } 

    /**
     * Ancestral conversion graph this list belongs to.
     */
    private final ConversionGraph acg;

    /**
     * List of events on clonal frame.
     */
    private final List<Event> events;
    private boolean dirty;

    public CFEventList(ConversionGraph acg) {
        this.acg = acg;
        
        events = new ArrayList<>();
        dirty = true;
    }

    /**
     * Obtain ordered list of events that make up the clonal frame.  Used
     * for ARG probability density calculations and for various state proposal
     * operators.
     * 
     * @return List of events.
     */
    public List<Event> getCFEvents() {
        updateEvents();
        
        return events;
    }

    /**
     * Mark the event list as dirty.
     */
    public void makeDirty() {
        dirty = true;
    }

    /**
     * Assemble sorted list of events on clonal frame and a map from nodes
     * to these events.
     */
    public void updateEvents() {
        if (!dirty)
            return;
        
        events.clear();
        
        // Create event list
        for (Node node : acg.getNodesAsArray()) {
            Event event = new Event(node);
            events.add(event);
        }
        
        // Sort events in increasing order of their times
        Collections.sort(events, (Event o1, Event o2) -> {
            if (o1.t<o2.t)
                return -1;
            
            if (o2.t<o1.t)
                return 1;
            
            return 0;
        });
        
        // Compute lineage counts:
        int k=0;
        for (Event event : events) {
            if (event.type == EventType.SAMPLE)
                k += 1;
            else
                k -= 1;
            
            event.lineages = k;
        }


        dirty = false;
    }
}
