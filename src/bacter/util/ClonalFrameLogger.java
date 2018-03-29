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

package bacter.util;

import bacter.ConversionGraph;
import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.io.PrintStream;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Logs clonal frame corresponding to conversion graph.")
public class ClonalFrameLogger extends CalculationNode implements Loggable {

    public Input<ConversionGraph> acgInput = new Input<>(
            "acg", "Conversion graph whose clonal frame you want to log.",
            Validate.REQUIRED);

    ConversionGraph acg;

    @Override
    public void initAndValidate() {
        acg = acgInput.get();
    }
    
    public void init(PrintStream out) {
        Node node = acg.getRoot();
        out.println("#NEXUS\n");
        out.println("Begin taxa;");
        out.println("\tDimensions ntax=" + acg.getLeafNodeCount() + ";");
        out.println("\t\tTaxlabels");
        acg.printTaxa(node, out, acg.getNodeCount() / 2);
        out.println("\t\t\t;");
        out.println("End;");

        out.println("Begin trees;");
        out.println("\tTranslate");
        acg.printTranslate(node, out, acg.getNodeCount() / 2);
        out.print(";");
    }

    @Override
    public void log(long sample, PrintStream out) {
        Tree tree = (Tree) acg.getCurrent();
        out.print("tree STATE_" + sample + " = ");
        // Don't sort, this can confuse CalculationNodes relying on the tree
        //tree.getRoot().sort();
        final int[] dummy = new int[1];
        final String newick = tree.getRoot().toSortedNewick(dummy, false);
        out.print(newick);
        out.print(";");
    }


    @Override
    public void close(PrintStream out) {
        out.print("End;");
    }

}
