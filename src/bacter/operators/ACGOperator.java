/*
 * Copyright (C) 2013 Tim Vaughan <tgvaughan@gmail.com>
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

package bacter.operators;

import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.Region;
import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Node;
import feast.input.In;

/**
 * Abstract class of operators which act on the ConversionGraph state.
 * 
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public abstract class ACGOperator extends Operator {
    
    public Input<ConversionGraph> argInput = new In<ConversionGraph>(
            "acg", "Ancestral conversion graph.").setRequired();
    
    protected ConversionGraph acg;

    protected class ProposalFailed extends Exception {
        public ProposalFailed() { }
    };

    @Override
    public void initAndValidate() throws Exception {
        acg = argInput.get();
    }
    
    /**
     * Return sister of node.
     * 
     * @param node
     * @return sister node
     */
    protected Node getSibling(Node node) {
        Node parent = node.getParent();
        if (parent.getLeft() == node)
            return parent.getRight();
        else
            return parent.getLeft();
    }

    /**
     * @return true if region map is valid under the restricted CO model.
     */
    protected boolean regionMapIsValid() {
        for (int c=0; c<acg.getConversions().size()-1; c++) {
            Conversion thisConv = acg.getConversions().get(c);
            Conversion nextConv = acg.getConversions().get(c+1);
            if (nextConv.getStartSite()-thisConv.getEndSite()<2)
                return false;
        }
        return true;
    }
}
