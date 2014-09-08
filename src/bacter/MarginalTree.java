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
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

/**
 * Light-weight Tree object representing marginal tree.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class MarginalTree {

    Node marginalRoot;

    public enum EventType { CF_LEAF, CF_COALESCENCE, CONV_DEPART, CONV_ARRIVE };
    private static class ACGEvent {
        public double time;
        public EventType type;
        public Node node;
        public Conversion conversion;
    }

    List<ACGEvent> eventList;
    
    public Node getRoot() {
        return marginalRoot;
    }

	public MarginalTree(Node marginalRoot) {
		this.marginalRoot = marginalRoot;
	}

    public MarginalTree(ConversionGraph acg, Region region) {
        this(acg, region.activeConversions);
    }

    public MarginalTree(ConversionGraph acg, Set<Conversion> convSet) {
        Set<Node> activeCFlineages = Sets.newHashSet();
        List<Node> inactiveCFlineages = Lists.newArrayList(acg.getExternalNodes());

        Set<Conversion> activeConversionLineages = Sets.newHashSet();

        inactiveCFlineages.sort((Node o1, Node o2) -> {
            if (o1.getHeight()<o2.getHeight())
                return -1;
            if (o1.getHeight()>o2.getHeight())
                return 1;
            return 0;
        });

        while (!inactiveCFlineages.isEmpty() &&
            activeCFlineages.size() + activeConversionLineages.size() > 1) {

            Conversion nextConversionStart = null;
            for (Conversion conv : convSet) {
                if (activeCFlineages.contains(conv.node1) &&
                    (nextConversionStart == null || conv.height1 < nextConversionStart.height1)) {
                    nextConversionStart = conv;
                }
            }

            Conversion nextConversionEnd = null;
            for (Conversion conv : convSet) {
                if (nextConversionEnd == null || conv.height2 < nextConversionEnd.height2) {
                    nextConversionEnd = conv;
                }
            }

            Node nextCFCoalescentNode = null;
            for (Node node : activeCFlineages) {
                if (nextCFCoalescentNode == null
                    || (!node.isRoot()
                    && node.getParent().getHeight() < nextCFCoalescentNode.getParent().getHeight())) {
                    nextCFCoalescentNode = node;
                }
            }

            Node nextActivationNode = null;
            if (!inactiveCFlineages.isEmpty())
                nextActivationNode = inactiveCFlineages.get(0);

        }
    }
}
