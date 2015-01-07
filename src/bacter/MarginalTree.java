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
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Light-weight Tree object representing marginal tree.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class MarginalTree {

    Node marginalRoot;

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

        Map<Conversion, Node> activeConversions = new HashMap<>();
        Map<Node, Node> activeCFlineages = new HashMap<>();

        for (ACGEventList.Event event : acg.getACGEvents()) {

            switch (event.getType()) {
                case CF_LEAF:
                    Node marginalLeaf = new Node();
                    marginalLeaf.setHeight(event.getTime());
                    marginalLeaf.setID(event.getNode().getID());
                    marginalLeaf.setNr(event.getNode().getNr());
                    activeCFlineages.put(event.getNode(), marginalLeaf);
                    break;

                case CF_COALESCENCE:
                    if (activeCFlineages.containsKey(event.getNode().getLeft())
                        && activeCFlineages.containsKey(event.getNode().getRight())) {
                        Node marginalNode = new Node();
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
                            Node marginalNode = activeCFlineages.get(event.getNode().getLeft());
                            activeCFlineages.remove(event.getNode().getLeft());
                            activeCFlineages.put(event.getNode(), marginalNode);
                        }

                        if (activeCFlineages.containsKey(event.getNode().getRight())) {
                            Node marginalNode = activeCFlineages.get(event.getNode().getRight());
                            activeCFlineages.remove(event.getNode().getRight());
                            activeCFlineages.put(event.getNode(), marginalNode);
                        }
                    }
                    break;

                case CONV_DEPART:
                    break;

                case CONV_ARRIVE:
                    break;
                default:
                    // Null type
                    break;
            }
            
        }
    }
}
