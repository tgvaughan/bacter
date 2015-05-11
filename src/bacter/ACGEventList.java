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

import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Node;
import java.util.ArrayList;
import java.util.List;

/**
 * Maintains an ordered list of events that make up the Ancestral Conversion
 * Graph. Used in the assembly of marginal trees.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ACGEventList {

    public enum EventType { CF_LEAF, CF_COALESCENCE, CONV_DEPART, CONV_ARRIVE }
    public static class Event {
        EventType type;
        double t;
        Node node;
        Conversion conversion;

        public EventType getType() {
            return type;
        }

        public double getTime() {
            return t;
        }

        public Node getNode() {
            return node;
        }

        public Conversion getConversion() {
            return conversion;
        }

        @Override
        public String toString() {
            return "t: " + t
                + " type: " + type
                + " node:" + node
                + " conv: " + conversion;
        }
    }

    /**
     * Event list.
     */
    private final List<Event> events; 

    /**
     * Event list dirty flag.
     */
    private boolean dirty;

    /**
     * Ancestral conversion graph this list belongs to.
     */
    private final ConversionGraph acg;

    /**
     * Construct a new event list for the given ACG.  There should only
     * be one of these objects per ACG object, created during the ACG
     * initAndValidate().
     * 
     * @param acg Conversion graph from which to compute event list.
     */
    public ACGEventList(ConversionGraph acg) {
        this.acg = acg;
        events = new ArrayList<>();
        dirty = true;
    }

    /**
     * Obtain sorted list of events that make up the ACG.
     * 
     * @return ACG event list.
     */
    public synchronized List<Event> getACGEvents() {
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
     * Assemble sorted list of events on ACG and a map from nodes and
     * conversions to these events.
     */
    public void updateEvents() {
        if (!dirty) {
            return;
        }

        events.clear();

        // Create unsorted event list.

        // Add CF events:
        for (Node node : acg.getNodesAsArray()) {
            Event event = new Event();
            event.node = node;
            if (node.isLeaf())
                event.type = EventType.CF_LEAF;
            else
                event.type = EventType.CF_COALESCENCE;
            event.t = node.getHeight();

            events.add(event);
        }

        // Add conversion events:
        for (Alignment alignment : acg.getAlignments()) {
            for (Conversion conv : acg.getConversions(alignment)) {
                Event departureEvent = new Event();
                departureEvent.conversion = conv;
                departureEvent.type = EventType.CONV_DEPART;
                departureEvent.t = conv.height1;
                events.add(departureEvent);

                Event arrivalEvent = new Event();
                arrivalEvent.conversion = conv;
                arrivalEvent.type = EventType.CONV_ARRIVE;
                arrivalEvent.t = conv.height2;
                events.add(arrivalEvent);
            }
        }

        // Sort event list:

        events.sort((Event o1, Event o2) -> {
            if (o1.t < o2.t)
                return -1;
            if (o1.t > o2.t)
                return 1;
            return 0;
        });

        dirty = false;
    }
    
}
