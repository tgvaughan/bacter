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

package bacter.util;

import bacter.ConversionGraph;
import bacter.Locus;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Class representing ACG log files.  Includes methods for
 * querying the number of ACGs defined, included and excluded
 * by the given burn-in percentage, as well as implementing an
 * iterator over all ACGs included after burn-in.  The iterator
 * automatically displays a progress bar on stdout.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class BacterACGLogReader implements ACGLogReader {
    File logFile;
    BufferedReader reader;

    List<String> preamble, postamble;
    String nextLine;

    List<Locus> loci;

    int nACGs, burnin;

    /**
     * Construct and initialize the reader.  The Preamble is
     * read and the list of loci constructed immediately.
     *
     * @param logFile ACG log file.
     * @throws IOException
     */
    public BacterACGLogReader(File logFile, double burninPercentage) throws IOException {
        this.logFile = logFile;

        reader = new BufferedReader(new FileReader(logFile));

        preamble = new ArrayList<>();
        skipPreamble();

        nACGs = 0;
        while (true) {
            if (getNextTreeString() == null)
                break;

            nACGs += 1;
        }
        burnin = (int)Math.round(nACGs*burninPercentage/100);

        postamble = new ArrayList<>();
        readPostamble();

        loci = new ArrayList<>();
        extractLoci();
    }


    /**
     * Internal method for skimming the preamble at the start
     * of the log, before we get to the tree section.
     *
     * @throws IOException
     */
    private void skipPreamble() throws IOException {
        boolean recordPreamble = preamble.isEmpty();

        while(true) {
            nextLine = reader.readLine();

            if (nextLine == null)
                throw new IOException("Reached end of file while searching for first tree.");

            nextLine = nextLine.trim();

            if (nextLine.toLowerCase().startsWith("tree"))
                break;

            if (recordPreamble)
                preamble.add(nextLine);
        }
    }

    /**
     * Internal method for extracting postamble following trees.
     *
     * @throws IOException
     */
    private void readPostamble() throws IOException {
        while (true) {
            if (nextLine == null)
                break;

            postamble.add(nextLine);

            nextLine = reader.readLine();
        }
    }

    /**
     * @return Everything read from the log file up until the first tree line.
     */
    public String getPreamble() {
        StringBuilder sb = new StringBuilder();
        for (String line : preamble)
            sb.append(line).append("\n");

        return sb.toString();
    }

    /**
     * @return Everything read from the log file following the last tree line.
     */
    public String getPostamble() {
        StringBuilder sb = new StringBuilder();
        for (String line : postamble)
            sb.append(line).append("\n");

        return sb.toString();
    }

    /**
     * Retrieve list of loci from preamble or postamble.
     */
    private void extractLoci() {
        List<String> prepost = new ArrayList<>();
        prepost.addAll(preamble);
        prepost.addAll(postamble);

        for (String line : prepost) {
            line = line.trim();
            if (line.startsWith("loci ") && line.endsWith(";")) {
                for (String locusEntry : line.substring(5,line.length()-1).split(" ")) {
                    String[] locusPair = locusEntry.split(":");
                    loci.add(new Locus(locusPair[0], Integer.parseInt(locusPair[1])));
                }
            }
        }
    }

    /**
     * Rewind to the beginning of the file.
     *
     * @throws IOException
     */
    private void reset() throws IOException {
        reader.close();
        reader = new BufferedReader(new FileReader(logFile));
        skipPreamble();
    }

    /**
     * @return the next available tree string or null if none exists
     * @throws IOException
     */
    private String getNextTreeString() throws IOException {
        StringBuilder sb = new StringBuilder();

        while (true) {
            if (nextLine == null || nextLine.trim().toLowerCase().equals("end;"))
                return null;

            sb.append(nextLine.trim());
            if (nextLine.trim().endsWith(";"))
                break;

            nextLine = reader.readLine();
        }
        nextLine = reader.readLine();

        String treeString = sb.toString();

        return treeString.substring(treeString.indexOf("("));
    }

    /**
     * Skip burn-in portion of log.
     *
     * @throws IOException
     */
    private void skipBurnin() throws IOException {
        for (int i=0; i<burnin; i++)
            getNextTreeString();
    }

    /**
     * @return loci read from the preamble
     */
    public List<Locus> getLoci() {
        return loci;
    }

    /**
     * @return total number of ACGs defined by file.
     */
    public int getACGCount() {
        return nACGs;
    }

    /**
     * @return number of ACGs excluded as burn-in
     */
    public int getBurnin() {
        return burnin;
    }

    /**
     * @return number of ACGs excluding burn-in
     */
    public int getCorrectedACGCount() {
        return nACGs - burnin;
    }

    /**
     * Retrieve an iterator for iterating over the ACGs represented
     * by this log file.  Important points
     *
     * 1. The iterator only iterates over as many (non-burnin) ACGs as exist
     * in the file when the ACGLogFileReader is constructed.  This is to avoid
     * problems associated with summarising ongoing analyses.
     *
     * 2. The iterator reuses a single ConversionGraph object during
     * the iteration.  This means that if you want to collect these
     * graphs as the iteration progresses you'll need to use
     * ConversionGraph::copy.
     *
     * @return ConversionGraph iterator
     */
    @Override
    public Iterator<ConversionGraph> iterator() {
        try {
            reset();
            skipBurnin();
        } catch (IOException e) {
            throw new IllegalStateException(e.getMessage());
        }

        ConversionGraph acg = new ConversionGraph();
        for (Locus locus : getLoci())
            acg.lociInput.setValue(locus, acg);
        try {
            acg.initAndValidate();
        } catch (Exception e) {
            throw new IllegalStateException(e.getMessage());
        }

        return new Iterator<ConversionGraph>() {

            boolean lineConsumed = true;
            String nextLine = null;

            int current = 0;

            private String getNextLineNoConsume() {
                if (lineConsumed) {
                    try {
                        nextLine = getNextTreeString();
                        lineConsumed = false;
                    } catch (IOException e) {
                        throw new IllegalStateException(e.getMessage());
                    }
                }

                return nextLine;
            }

            private void printProgressBar() {

                if (current==0) {
                    System.out.println("0%             25%            50%            75%           100%");
                    System.out.println("|--------------|--------------|--------------|--------------|");
                }

                if (current < getCorrectedACGCount()-1) {
                    if (current % (int) Math.ceil(getCorrectedACGCount() / 61.0) == 0) {
                        System.out.print("\r");
                        for (int i = 0; i < Math.round(61.0 * current / getCorrectedACGCount()); i++)
                            System.out.print("*");
                        System.out.flush();
                    }
                } else {
                    System.out.print("\r");
                    for (int i=0; i<61; i++)
                        System.out.print("*");
                    System.out.println();
                }

            }

            @Override
            public boolean hasNext() {
                return current<getCorrectedACGCount() && getNextLineNoConsume() != null;
            }

            @Override
            public ConversionGraph next() {
                String result = getNextLineNoConsume();
                lineConsumed = true;
                acg.fromExtendedNewick(result);

                printProgressBar();
                current += 1;

                return acg;
            }
        };
    }
}

