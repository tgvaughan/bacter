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
import bacter.util.BacterACGLogReader;
import bacter.util.COACGLogFileReader;
import beast.evolution.tree.Node;

import javax.xml.stream.XMLStreamException;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class DifferenceFromTrueACG {

    private static class Options {
        double boundaryTol = 0.2;
        double ageTol = 0.2;
        File logFile, truthFile, outFile;
        boolean useCOFormat = false;
    }

    public static void printUsageAndExit(int exitCode) {
        System.out.println("Usage: DifferenceFromTrueACG [-boundaryTol t] [-ageTol t] [-co] truth.tree log.trees output_file");
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
                case "co":
                    options.useCOFormat = true;
                    break;

                case "boundaryTol":
                    i += 1;
                    if (i>=args.length)
                        printUsageAndExit(1);
                    try {
                        options.boundaryTol = Double.valueOf(args[i]);
                    } catch (NumberFormatException e) {
                        System.out.println("Argument to -boundaryTol must be a number.");
                        printUsageAndExit(1);
                    }
                    break;

                case "ageTol":
                    i += 1;
                    if (i>=args.length)
                        printUsageAndExit(1);
                    try {
                        options.ageTol = Double.valueOf(args[i]);
                    } catch (NumberFormatException e) {
                        System.out.println("Argument to -ageTol must be a number.");
                        printUsageAndExit(1);
                    }
                    break;

                default:
                    System.err.println("Unknown argument: " + args[i]);
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

    public static class Clade extends BitSet {
        public double age;
    };

    public static Clade getClades(Clade[] clades, Node node) {
        Clade clade = new Clade();
        clade.age = node.getHeight();

        if (node.isLeaf()) {
            clade.set(node.getNr());
        } else {
            for (Node child : node.getChildren())
                clade.or(getClades(clades, child));
        }

        clades[node.getNr()] = clade;

        return clade;
    }

    /**
     * Count number of true clades which exist in the provided sampled ARG.
     *
     * @param trueClades clades in true arg
     * @param clades clades in sampled arg
     * @param ageTol maximum relative age error
     * @param cladeHist
     * @return number of found clades
     */
    public static int countFoundClades(Clade[] trueClades, Clade[] clades, double ageTol,
                                       Map<Clade, Integer> cladeHist) {
        int foundClades = 0;
        for (Clade trueClade : trueClades) {
            for (Clade clade : clades) {
                if (!clade.equals(trueClade) || (clade.cardinality()>1 && (trueClade.age - clade.age)/trueClade.age > ageTol))
                    continue;

                foundClades += 1;

                cladeHist.put(clade, cladeHist.get(clade)+1);

                break;
            }
        }

        return foundClades;
    }

    /**
     * Count the number of true conversions which have correspondences on the
     * provided sampled ARG.
     *
     * @param trueACG true arg
     * @param trueClades clades in true arg
     * @param acg sampled arg
     * @param clades clades in sampled arg
     * @param boundaryTol minimum relative error in region boundaries to allow.
     * @param convHist
     * @return number of found conversions
     */
    public static int countFoundConversions(ConversionGraph trueACG, Clade[] trueClades,
                                            ConversionGraph acg, Clade[] clades,
                                            double boundaryTol, double ageTol,
                                            Map<Conversion, Integer> convHist) {
        int count = 0;
        for (Conversion trueConv : trueACG.getConversions(trueACG.getLoci().get(0))) {
            Clade trueFromClade = trueClades[trueConv.getNode1().getNr()];
            Clade trueToClade = trueClades[trueConv.getNode1().getNr()];
            for (Conversion conv : acg.getConversions(acg.getLoci().get(0)))  {
                Clade fromClade = clades[conv.getNode1().getNr()];
                Clade toClade = clades[conv.getNode1().getNr()];

                if (fromClade.equals(trueFromClade) && toClade.equals(trueToClade)) {
                    if (    Math.abs(trueConv.getStartSite()-conv.getStartSite())/trueConv.getSiteCount() <= boundaryTol &&
                            Math.abs(trueConv.getEndSite()-conv.getEndSite())/trueConv.getSiteCount() <= boundaryTol &&
//                            Math.abs(trueConv.getHeight1()-conv.getHeight1())/trueConv.getHeight1() <= ageTol &&
                            Math.abs(trueConv.getHeight2()-conv.getHeight2())/trueConv.getHeight2() <= ageTol)
                    {
                        count += 1;

                        convHist.put(trueConv, convHist.get(trueConv) + 1);

                        break;
                    }
                }
            }
        }

        return count;
    }

    public static void main(String[] args) throws IOException, XMLStreamException {

        Options options = processArguments(args);

        // Load true ARG

        BacterACGLogReader truthReader = new BacterACGLogReader(options.truthFile, 0);
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

         // Set up histograms

        Map<Clade, Integer> cladeHist = new HashMap<>();
        for (Clade clade : trueClades)
            cladeHist.put(clade, 0);

        Map<Conversion, Integer> convHist = new HashMap<>();
        for (Conversion conv : trueACG.getConversions(trueACG.getLoci().get(0)))
            convHist.put(conv, 0);

        // Set up ARG log file reader

        Iterable<ConversionGraph> logReader;
        if (options.useCOFormat) {
            logReader = new COACGLogFileReader(options.logFile, 0);
        } else {
            logReader = new BacterACGLogReader(options.logFile, 0);
        }


        // Compute and write summary statistics to output file

        try (PrintStream ps = new PrintStream(options.outFile)) {
            ps.println("trueCladeCount recoveredCladeCount trueConvCount sampledConvCount recoveredConvCount meanTimeError maxTimeError");

            for (ConversionGraph acg : logReader) {

                Clade[] clades = new Clade[acg.getNodeCount()];
                getClades(clades, acg.getRoot());

                List<Double> timeErrors = new ArrayList<>();
                int foundClades = countFoundClades(trueClades, clades, options.ageTol, cladeHist);
                int foundConvs = countFoundConversions(trueACG, trueClades, acg, clades,
                        options.boundaryTol, options.ageTol, convHist);

                ps.print(trueACG.getNodeCount() + "\t" +
                        foundClades + "\t" +
                        trueACG.getConvCount(trueACG.getLoci().get(0)) + "\t" +
                        acg.getConvCount(acg.getLoci().get(0)) + "\t" +
                        foundConvs + "\n");
            }
        }



    }
}
