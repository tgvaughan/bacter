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
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

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

        for (ACGEventList.Event event : acg.getACGEvents()) {

            switch (event.getType()) {
                case CF_LEAF:
                    MarginalNode marginalLeaf = new MarginalNode();
                    marginalLeaf.setHeight(event.getTime());
                    marginalLeaf.setID(event.getNode().getID());
                    marginalLeaf.setNr(event.getNode().getNr());
                    marginalLeaf.cfNodeNr = event.getNode().getNr();
                    activeCFlineages.put(event.getNode(), marginalLeaf);
                    break;

                case CF_COALESCENCE:
                    if (activeCFlineages.containsKey(event.getNode().getLeft())
                        && activeCFlineages.containsKey(event.getNode().getRight())) {

                        MarginalNode marginalNode = new MarginalNode();
                        marginalNode.setNr(nextNonLeafNr++);
                        marginalNode.cfNodeNr = event.getNode().getNr();
                        Node marginalLeft = activeCFlineages.get(event.getNode().getLeft());
                        Node marginalRight = activeCFlineages.get(event.getNode().getRight());

                        marginalNode.setHeight(event.getTime());
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

                case CONV_DEPART:
                    if (!convSet.contains(event.getConversion()))
                        break;

                    if (activeCFlineages.containsKey(event.getConversion().getNode1())) {
                        MarginalNode marginalNode = activeCFlineages.get(event.getConversion().getNode1());
                        activeCFlineages.remove(event.getConversion().getNode1());
                        activeConversions.put(event.getConversion(), marginalNode);
                    }
                    break;

                case CONV_ARRIVE:
                    if (!convSet.contains(event.getConversion()))
                        break;

                    if (activeCFlineages.containsKey(event.getConversion().getNode2())

                        && activeConversions.containsKey(event.getConversion())) {
                        MarginalNode marginalNode = new MarginalNode();
                        marginalNode.setNr(nextNonLeafNr++);
                        MarginalNode marginalLeft = activeCFlineages.get(event.getConversion().getNode2());
                        MarginalNode marginalRight = activeConversions.get(event.getConversion());

                        marginalNode.setHeight(event.getTime());
                        marginalNode.addChild(marginalLeft);
                        marginalNode.addChild(marginalRight);

                        activeConversions.remove(event.getConversion());
                        activeCFlineages.put(event.getConversion().getNode2(), marginalNode);

                    } else {

                        if (activeConversions.containsKey(event.getConversion())) {
                            MarginalNode marginalNode = activeConversions.get(event.getConversion());
                            activeConversions.remove(event.getConversion());
                            activeCFlineages.put(event.getConversion().getNode2(), marginalNode);
                        }
                    }
                    break;
                default:
                    // Null type
                    break;
            }

            // A single active CF lineage should remain:
            marginalRoot = activeCFlineages.get(acg.getRoot());
        }
    }

    @Override
    public String toString() {
        return marginalRoot.toString();
    }
}
