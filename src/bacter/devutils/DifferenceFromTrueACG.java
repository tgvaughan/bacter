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

import bacter.ConversionGraph;
import bacter.acgannotator.ACGLogFileReader;
import beast.app.util.Arguments;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class DifferenceFromTrueACG {

    private static class Options {
        double burnin = 10.0;
        File logFile, truthFile, outFile = new File("out.log");
        boolean useCOFormat = false;
    }

    public static void printUsageAndExit(int exitCode) {
        System.out.println("Usage: DifferenceFromTrueACG [-burnin n] [-co] truth.tree log.trees");
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

                default:
                    printUsageAndExit(1);
            }

            i++;
        }

        if (args.length-i < 2 || args.length-i > 3)
            printUsageAndExit(0);

        options.truthFile = new File(args[i++]);
        options.logFile = new File(args[i++]);

        if (i<args.length)
            options.outFile = new File(args[i]);

        return options;
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

        // Load ARGs from log file

        List<ConversionGraph> sampledACGs = new ArrayList<>();
        if (options.useCOFormat) {

        } else {
            ACGLogFileReader logReader = new ACGLogFileReader(options.logFile, options.burnin);

            for (ConversionGraph acg : logReader)
                sampledACGs.add(acg.copy());
        }

        // Compute difference summaries and write to output

    }
}
