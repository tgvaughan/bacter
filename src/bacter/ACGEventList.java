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
        public final EventType type;
        public final double t;
        public final Node node;
        public final Conversion conversion;

        public Event(Node cfNode) {
            node = cfNode;
            t = cfNode.getHeight();
            if (node.isLeaf())
                type = EventType.CF_LEAF;
            else
                type = EventType.CF_COALESCENCE;

            conversion = null;
        }

        public Event(Conversion conv, boolean isDeparture) {
            conversion = conv;
            if (isDeparture) {
                type = EventType.CONV_DEPART;
                node = conv.getNode1();
                t = conv.getHeight1();
            } else {
                type = EventType.CONV_ARRIVE;
                node = conv.getNode2();
                t = conv.getHeight2();
            }
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
     * Construct a new event list for the given ACG.  There should only
     * be one of these objects per ACG object, created during the ACG
     * initAndValidate().
     * 
     * @param acg Conversion graph from which to compute event list.
     */
    public ACGEventList(ConversionGraph acg) {
        this(acg, null);
    }

    /**
     * Construct a new event list for the given ACG and locus.
     *
     * @param acg Conversion graph from which to compute event list
     * @param locus conversions for this event only are included
     */
    public ACGEventList(ConversionGraph acg, Locus locus) {
        this.events = new ArrayList<>();

        // Create unsorted event list.

        // Add CF events:
        for (Node node : acg.getNodesAsArray()) {
            events.add(new Event(node));
        }

        // Add conversion events:
        if (locus == null) {
            for (Locus l : acg.getConvertibleLoci()) {
                for (Conversion conv : acg.getConversions(l)) {
                    events.add(new Event(conv, true));
                    events.add(new Event(conv, false));
                }
            }
        } else {
            for (Conversion conv : acg.getConversions(locus)) {
                events.add(new Event(conv, true));
                events.add(new Event(conv, false));
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
    }

    /**
     * Obtain sorted list of events that make up the ACG.
     * 
     * @return ACG event list.
     */
    public synchronized List<Event> getACGEvents() {
        return events;
    }
}
