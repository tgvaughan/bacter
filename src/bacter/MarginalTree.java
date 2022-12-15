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

import java.util.*;

/**
 * Light-weight Tree object representing marginal tree.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class MarginalTree {

    MarginalNode marginalRoot;

    public MarginalNode getRoot() {
        return marginalRoot;
    }

    public MarginalTree(ConversionGraph acg, Region region) {
        this(acg, region.activeConversions);
    }

    public MarginalTree(ConversionGraph acg, Set<Conversion> convSet) {

        Map<Conversion, MarginalNode> activeConversions = new HashMap<>();
        Map<Node, MarginalNode> activeCFlineages = new HashMap<>();

        int nextNonLeafNr = acg.getLeafNodeCount();

        ArrayList<ConversionEvent> convEvents = new ArrayList<>();
        convSet.forEach(conversion -> {
            convEvents.add(new ConversionEvent(conversion, false));
            convEvents.add(new ConversionEvent(conversion, true));
        });
        convEvents.sort((c1, c2) -> {
            if (c1.height < c2.height)
                return -1;

            if (c1.height > c2.height)
                return 1;

            return 0;
        });

        List<CFEventList.Event> cfEvents = acg.getCFEvents();

        int convEventIdx=0;
        for (int eventIdx=0; eventIdx<cfEvents.size(); eventIdx++) {
            CFEventList.Event event = cfEvents.get(eventIdx);

            switch (event.getType()) {
                case SAMPLE:
                    MarginalNode marginalLeaf = new MarginalNode();
                    marginalLeaf.setHeight(event.getHeight());
                    marginalLeaf.setID(event.getNode().getID());
                    marginalLeaf.setNr(event.getNode().getNr());
                    marginalLeaf.cfNodeNr = event.getNode().getNr();
                    activeCFlineages.put(event.getNode(), marginalLeaf);
                    break;

                case COALESCENCE:
                    if (activeCFlineages.containsKey(event.getNode().getLeft())
                            && activeCFlineages.containsKey(event.getNode().getRight())) {

                        MarginalNode marginalNode = new MarginalNode();
                        marginalNode.setNr(nextNonLeafNr++);
                        marginalNode.cfNodeNr = event.getNode().getNr();
                        Node marginalLeft = activeCFlineages.get(event.getNode().getLeft());
                        Node marginalRight = activeCFlineages.get(event.getNode().getRight());

                        marginalNode.setHeight(event.getHeight());
                        marginalNode.addChild(marginalLeft);
                        marginalNode.addChild(marginalRight);

                        activeCFlineages.remove(event.getNode().getLeft());
                        activeCFlineages.remove(event.getNode().getRight());
                        activeCFlineages.put(event.getNode(), marginalNode);

                    } else {

                        if (activeCFlineages.containsKey(event.getNode().getLeft())) {
                            MarginalNode marginalNode = activeCFlineages.get(event.getNode().getLeft());
                            activeCFlineages.remove(event.getNode().getLeft());
                            activeCFlineages.put(event.getNode(), marginalNode);
                            break;
                        }

                        if (activeCFlineages.containsKey(event.getNode().getRight())) {
                            MarginalNode marginalNode = activeCFlineages.get(event.getNode().getRight());
                            activeCFlineages.remove(event.getNode().getRight());
                            activeCFlineages.put(event.getNode(), marginalNode);
                            break;
                        }
                    }
                    break;
            }

            while (convEventIdx < convEvents.size() &&
                    (event.node.isRoot() || convEvents.get(convEventIdx).height < cfEvents.get(eventIdx + 1).getHeight())) {
                ConversionEvent convEvent = convEvents.get(convEventIdx++);

                if (convEvent.isDeparture) {
                    if (activeCFlineages.containsKey(convEvent.conversion.getNode1())) {
                        MarginalNode marginalNode = activeCFlineages.get(convEvent.conversion.getNode1());
                        activeCFlineages.remove(convEvent.conversion.getNode1());
                        activeConversions.put(convEvent.conversion, marginalNode);
                    }

                } else {
                    if (activeCFlineages.containsKey(convEvent.conversion.getNode2())

                            && activeConversions.containsKey(convEvent.conversion)) {
                        MarginalNode marginalNode = new MarginalNode();
                        marginalNode.setNr(nextNonLeafNr++);
                        MarginalNode marginalLeft = activeCFlineages.get(convEvent.conversion.getNode2());
                        MarginalNode marginalRight = activeConversions.get(convEvent.conversion);

                        marginalNode.setHeight(convEvent.height);
                        marginalNode.addChild(marginalLeft);
                        marginalNode.addChild(marginalRight);

                        activeConversions.remove(convEvent.conversion);
                        activeCFlineages.put(convEvent.conversion.getNode2(), marginalNode);

                    } else {

                        if (activeConversions.containsKey(convEvent.conversion)) {
                            MarginalNode marginalNode = activeConversions.get(convEvent.conversion);
                            activeConversions.remove(convEvent.conversion);
                            activeCFlineages.put(convEvent.conversion.getNode2(), marginalNode);
                        }
                    }
                }
            }
        }

        // A single active CF lineage should remain:
        marginalRoot = activeCFlineages.get(acg.getRoot());
    }

    @Override
    public String toString() {
        return marginalRoot.toString();
    }
}
