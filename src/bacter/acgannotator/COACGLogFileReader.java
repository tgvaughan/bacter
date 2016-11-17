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

package bacter.acgannotator;

import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.Locus;
import beast.evolution.tree.Node;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class COACGLogFileReader implements Iterable<ConversionGraph> {
    File logFile;

    XMLStreamReader xmlStreamReader;

    List<Locus> loci = new ArrayList<>();

    int nACGs, burnin;

    public COACGLogFileReader(File logFile, double burninPercentage) throws IOException, XMLStreamException {
        this.logFile = logFile;

        reset();

        List<Integer> locusSize = null;
        List<String> locusName = null;

        String text;

        nACGs = 0;

        while (xmlStreamReader.hasNext()) {
            int eventType = xmlStreamReader.next();
            if (eventType == XMLStreamConstants.START_ELEMENT) {
                switch (xmlStreamReader.getLocalName().toLowerCase()) {
                    case "iteration":
                        nACGs += 1;
                        break;

                    case "blocks":
                        xmlStreamReader.next();
                        text = xmlStreamReader.getText().trim();
                        locusSize = new ArrayList<>();
                        for (String rangeStr : text.split(";")) {
                            String[] splitRange = rangeStr.split(",");
                            locusSize.add(Integer.valueOf(splitRange[1]) - Integer.valueOf(splitRange[0]));
                        }
                        break;

                    case "regions":
                        xmlStreamReader.next();
                        locusName = Arrays.asList(xmlStreamReader.getText().trim().split(","));
                        break;

                    default:
                }
            }
        }

        if (locusName == null || locusSize == null)
            throw new IOException("Missing <blocks> or <regions> elements.");

        for (int i=0; i<locusName.size(); i++)
            loci.add(new Locus(locusName.get(i), locusSize.get(i)));

        burnin = (int)Math.round(nACGs*burninPercentage/100);
    }

    private void reset() throws IOException {
        try {
            if (xmlStreamReader != null)
                xmlStreamReader.close();

            xmlStreamReader = XMLInputFactory.newFactory().createXMLStreamReader(
                    new BufferedInputStream(new FileInputStream(logFile)));
        } catch (XMLStreamException | FileNotFoundException e) {
            throw new IOException(e.getMessage());
        }
    }

    /**
     * Move through log file until the opening tag of the named element.
     *
     * @param elementName name of element to skip until
     * @return true if element found, false otherwise
     * @throws XMLStreamException
     */
    public boolean skipUntil(String elementName) throws XMLStreamException {
        while (xmlStreamReader.hasNext()) {
            if (xmlStreamReader.next() == XMLStreamConstants.START_ELEMENT &&
                    xmlStreamReader.getLocalName().toLowerCase().equals(elementName))
                return true;
        }

        return false;
    }

    /**
     * Skip burnin portion of log file
     *
     * @throws XMLStreamException
     */
    public void skipBurnin() throws XMLStreamException {
        for (int i=0; i<burnin; i++)
            skipUntil("iteration");
    }

    @Override
    public Iterator<ConversionGraph> iterator() {

        try {
            reset();
            skipBurnin();
        } catch (XMLStreamException | IOException e) {
            throw new IllegalStateException(e.getMessage());
        }

        return new Iterator<ConversionGraph>() {

            int current = 0;

            @Override
            public boolean hasNext() {
                return current < nACGs;
            }

            @Override
            public ConversionGraph next() {

                current += 1;

                ConversionGraph acg = new ConversionGraph();
                for (Locus locus : loci)
                    acg.lociInput.setValue(locus, acg);
                try {
                    acg.initAndValidate();
                } catch (Exception e) {
                    throw new IllegalStateException(e.getMessage());
                }

                String newick = null;
                List<Integer> rStarts = new ArrayList<>();
                List<Integer> rEnds = new ArrayList<>();
                List<Integer> rEFroms = new ArrayList<>();
                List<Integer> rETos = new ArrayList<>();
                List<Double> rAFroms = new ArrayList<>();
                List<Double> rATos = new ArrayList<>();

                try {
                    skipUntil("iteration");

                    while (xmlStreamReader.hasNext()) {
                        int type = xmlStreamReader.next();

                        if (type == XMLStreamConstants.END_ELEMENT &&
                                xmlStreamReader.getLocalName().toLowerCase().equals("iteration"))
                            break;

                        if (type != XMLStreamConstants.START_ELEMENT)
                            continue;


                        switch (xmlStreamReader.getLocalName().toLowerCase()) {
                            case "tree":
                                xmlStreamReader.next();
                                newick = xmlStreamReader.getText().trim();
                                break;

                            case "recedge":

                                while ((type = xmlStreamReader.next()) != XMLStreamConstants.END_ELEMENT ||
                                        !xmlStreamReader.getLocalName().toLowerCase().equals("recedge")) {
                                    if (type != XMLStreamConstants.START_ELEMENT)
                                        continue;

                                    switch (xmlStreamReader.getLocalName().toLowerCase()) {
                                        case "start":
                                            xmlStreamReader.next();
                                            rStarts.add(Integer.valueOf(xmlStreamReader.getText().trim()));
                                            break;

                                        case "end":
                                            xmlStreamReader.next();
                                            rEnds.add(Integer.valueOf(xmlStreamReader.getText().trim()) - 1);
                                            break;

                                        case "efrom":
                                            xmlStreamReader.next();
                                            rEFroms.add(Integer.valueOf(xmlStreamReader.getText().trim()));
                                            break;

                                        case "eto":
                                            xmlStreamReader.next();
                                            rETos.add(Integer.valueOf(xmlStreamReader.getText().trim()));
                                            break;

                                        case "afrom":
                                            xmlStreamReader.next();
                                            rAFroms.add(Double.valueOf(xmlStreamReader.getText().trim()));
                                            break;

                                        case "ato":
                                            xmlStreamReader.next();
                                            rATos.add(Double.valueOf(xmlStreamReader.getText().trim()));
                                            break;

                                        default:
                                    }
                                }

                                break;

                            default:
                        }

                    }
                } catch (XMLStreamException e) {
                    throw new IllegalStateException(e.getMessage());
                }

                acg.fromExtendedNewick(newick, true, 0);
                for (int i=0; i<rStarts.size(); i++) {
                    Node fromNode = acg.getNode(rEFroms.get(i));
                    Node toNode = acg.getNode(rETos.get(i));

                    Conversion conv = new Conversion(
                            toNode,
                            toNode.getHeight() + rATos.get(i),
                            fromNode,
                            fromNode.getHeight() + rAFroms.get(i),
                            rStarts.get(i),
                            rEnds.get(i),
                            acg,
                            loci.get(0));
                    acg.addConversion(conv);
                }

                current += 1;

                return acg;
            }
        };
    }

    public static void main(String[] args) throws IOException, XMLStreamException {

        COACGLogFileReader reader = new COACGLogFileReader(
                new File("/home/tvaughan/articles/bacter-paper/simulation_studies/robustness/blah.xml"),
                0);

        System.out.println("Iterations: " + reader.nACGs);

        for (ConversionGraph acg : reader)
            System.out.println(acg);
    }
}
