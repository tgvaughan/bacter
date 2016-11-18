/*
 * Copyright (C) 2016 Tim Vaughan <tgvaughan@gmail.com>
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

package bacter.devutils;

import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.acgannotator.ACGLogFileReader;
import beast.evolution.tree.Node;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class DifferenceFromTrueACG {

    private static class Options {
        double burnin = 10.0;
        double overlapTol = 50.0;
        File logFile, truthFile, outFile;
        boolean useCOFormat = false;
    }

    public static void printUsageAndExit(int exitCode) {
        System.out.println("Usage: DifferenceFromTrueACG [-burnin n] [-overlapTol] [-co] truth.tree log.trees output_file");
        System.exit(exitCode);
    }

    /**
     * Process command line arguments.
     *
     * @param args
     * @return
     */
    public static Options processArguments(String[] args) {

        Options options = new Options();

        int i=0;
        while (i<args.length && args[i].startsWith("-")) {
            switch (args[i].substring(1)) {
                case "burnin":
                    i += 1;
                    if (i>=args.length)
                        printUsageAndExit(1);
                    try {
                        options.burnin = Double.valueOf(args[i]);
                    } catch (NumberFormatException e) {
                        System.out.println("Argument to -burnin must be a number.");
                        printUsageAndExit(1);
                    }
                    break;

                case "co":
                    options.useCOFormat = true;
                    break;

                case "overlapTol":
                    i += 1;
                    if (i>=args.length)
                        printUsageAndExit(1);
                    try {
                        options.overlapTol = Double.valueOf(args[i]);
                    } catch (NumberFormatException e) {
                        System.out.println("Argument to -overlapTol must be a number.");
                        printUsageAndExit(1);
                    }

                default:
                    printUsageAndExit(1);
            }

            i++;
        }

        if (args.length-i < 3)
            printUsageAndExit(0);

        options.truthFile = new File(args[i++]);
        options.logFile = new File(args[i++]);
        options.outFile = new File(args[i]);

        return options;
    }

    public static class Clade extends BitSet { };
    public static class CladePair {
        Clade from, to;

        public CladePair(Clade from, Clade to) {
            this.from = from;
            this.to = to;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            CladePair cladePair = (CladePair) o;

            return from.equals(cladePair.from) && to.equals(cladePair.to);
        }

        @Override
        public int hashCode() {
            int result = from.hashCode();
            result = 31 * result + to.hashCode();
            return result;
        }
    }

    public static Clade getClades(Clade[] clades, Node node) {
        Clade clade = new Clade();

        if (node.isLeaf()) {
            clade.set(node.getNr());
        } else {
            for (Node child : node.getChildren())
                clade.or(getClades(clades, child));
        }

        clades[node.getNr()] = clade;

        return clade;
    }

    public static void main(String[] args) throws IOException {

        Options options = processArguments(args);

        // Load true ARG

        ACGLogFileReader truthReader = new ACGLogFileReader(options.truthFile, 0);
        if (truthReader.getACGCount() != 1) {
            System.out.println("Expected exactly 1 ACG in truth file. Found " +
                    truthReader.getACGCount());
            System.exit(1);
        }

        ConversionGraph trueACG = null;
        for (ConversionGraph acg : truthReader)
            trueACG = acg;

        // Determine clades present in truth
        Clade[] trueClades = new Clade[trueACG.getNodeCount()];
        getClades(trueClades, trueACG.getRoot());
        Set<Clade> trueCladeSet = new HashSet<>(Arrays.asList(trueClades));

        // Set up ARG log file reader

        ACGLogFileReader logReader;
        if (options.useCOFormat) {
            throw new IllegalStateException("Unimplemented!");
        } else {
            logReader = new ACGLogFileReader(options.logFile, options.burnin);
        }

        // Compute and write summary statistics to output file

        try (PrintStream ps = new PrintStream(options.outFile)) {
            ps.println("cladeCountError trueConvCount recoveredConvCount");

            for (ConversionGraph acg : logReader) {

                Clade[] clades = new Clade[acg.getNodeCount()];
                getClades(clades, acg.getRoot());
                Set<Clade> cladeSet = new HashSet<>(Arrays.asList(clades));

                int foundConvs = 0;
                for (Conversion trueConv : trueACG.getConversions(trueACG.getLoci().get(0))) {
                    Clade trueFromClade = trueClades[trueConv.getNode1().getNr()];
                    Clade trueToClade = trueClades[trueConv.getNode1().getNr()];
                    for (Conversion conv : acg.getConversions(acg.getLoci().get(0)))  {
                        Clade fromClade = clades[conv.getNode1().getNr()];
                        Clade toClade = clades[conv.getNode1().getNr()];

                        if (fromClade.equals(trueFromClade) && toClade.equals(trueToClade)) {
                            int overlap = Math.min(conv.getEndSite(), trueConv.getEndSite()) -
                                    Math.max(conv.getStartSite(), trueConv.getStartSite());

                            if (overlap/(double)trueConv.getSiteCount()>options.overlapTol/100.0) {
                                foundConvs += 1;
                                break;
                            }
                        }
                    }
                }


                ps.print(trueACG.getConvCount(trueACG.getLoci().get(0)) + "\t" +
                        acg.getConvCount(acg.getLoci().get(0)) + "\t" +
                        foundConvs + "\n");
            }
        }



    }
}
