/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
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

package bacter.xmltests;

import beast.util.Randomizer;
import beast.util.XMLParser;
import org.junit.Test;
import test.beast.beast2vs1.trace.Expectation;
import test.beast.beast2vs1.trace.LogAnalyser;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import static junit.framework.Assert.assertTrue;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SimulatedACGTest {

    @Test
    public void test2Taxon() throws Exception {
        Randomizer.setSeed(1);

        XMLParser parser = new XMLParser();
        beast.core.Runnable runnable = parser.parseFile(
                new File("examples/ACGsimulations/simulateACGs2taxon.xml"));
        runnable.run();

        List<Expectation> expectations = new ArrayList<>();
        expectations.add(new Expectation("acg.CFheight", 1.0, 1e-2));
        expectations.add(new Expectation("acg.CFlength", 2.0, 1e-2));
        expectations.add(new Expectation("acg.nConv", 10.0, 5e-2));

        LogAnalyser logAnalyser = new LogAnalyser("simulateACGs2taxon.stats",
                expectations);

        for (int i=0; i<expectations.size(); i++) {
            assertTrue(expectations.get(i).isValid());
            assertTrue(expectations.get(i).isPassed());
        }

        Files.deleteIfExists(Paths.get("simulateACGs2taxon.stats"));
        Files.deleteIfExists(Paths.get("simulateACGs2taxon.converted"));
        Files.deleteIfExists(Paths.get("simulateACGs2taxon.trees"));
    }

    @Test
    public void test5Taxon() throws Exception {
        Randomizer.setSeed(1);

        XMLParser parser = new XMLParser();
        beast.core.Runnable runnable = parser.parseFile(
                new File("examples/ACGsimulations/simulateACGs5taxon.xml"));
        runnable.run();

        List<Expectation> expectations = new ArrayList<>();
        expectations.add(new Expectation("acg.CFheight", 1.606, 1e-2));
        expectations.add(new Expectation("acg.CFlength", 4.181, 1e-2));
        expectations.add(new Expectation("acg.nConv", 21.0, 5e-2));

        LogAnalyser logAnalyser = new LogAnalyser("simulateACGs5taxon.stats",
                expectations);

        for (int i=0; i<expectations.size(); i++) {
            assertTrue(expectations.get(i).isValid());
            assertTrue(expectations.get(i).isPassed());
        }

        Files.deleteIfExists(Paths.get("simulateACGs5taxon.stats"));
        Files.deleteIfExists(Paths.get("simulateACGs5taxon.converted"));
        Files.deleteIfExists(Paths.get("simulateACGs5taxon.trees"));
    }

    @Test
    public void test5TaxonDynamicPopSize() throws Exception {
        Randomizer.setSeed(1);

        XMLParser parser = new XMLParser();
        beast.core.Runnable runnable = parser.parseFile(
                new File("examples/ACGsimulations/simulateACGs5taxonDynamicPopSize.xml"));
        runnable.run();

        List<Expectation> expectations = new ArrayList<>();
        expectations.add(new Expectation("acg.CFheight", 8.840, 1e-2));
        expectations.add(new Expectation("acg.CFlength", 25.312, 1e-2));
        expectations.add(new Expectation("acg.nConv", 25.464, 5e-2));

        LogAnalyser logAnalyser = new LogAnalyser("simulateACGs5taxonDynamicPopSize.stats",
                expectations);

        for (int i=0; i<expectations.size(); i++) {
            assertTrue(expectations.get(i).isValid());
            assertTrue(expectations.get(i).isPassed());
        }

        Files.deleteIfExists(Paths.get("simulateACGs5taxonDynamicPopSize.stats"));
        Files.deleteIfExists(Paths.get("simulateACGs5taxonDynamicPopSize.converted"));
        Files.deleteIfExists(Paths.get("simulateACGs5taxonDynamicPopSize.trees"));
    }

    @Test
    public void test5TaxonSerialSampling() throws Exception {
        Randomizer.setSeed(1);

        XMLParser parser = new XMLParser();
        beast.core.Runnable runnable = parser.parseFile(
                new File("examples/ACGsimulations/simulateACGs5taxonSerialSampling.xml"));
        runnable.run();

        List<Expectation> expectations = new ArrayList<>();
        expectations.add(new Expectation("acg.CFheight", 1.909, 1e-2));
        expectations.add(new Expectation("acg.CFlength", 4.655, 1e-2));
        expectations.add(new Expectation("acg.nConv", 23.381, 5e-2));

        LogAnalyser logAnalyser = new LogAnalyser("simulateACGs5taxonSerialSampling.stats",
                expectations);

        for (int i=0; i<expectations.size(); i++) {
            assertTrue(expectations.get(i).isValid());
            assertTrue(expectations.get(i).isPassed());
        }

        Files.deleteIfExists(Paths.get("simulateACGs5taxonSerialSampling.stats"));
        Files.deleteIfExists(Paths.get("simulateACGs5taxonSerialSampling.converted"));
        Files.deleteIfExists(Paths.get("simulateACGs5taxonSerialSampling.trees"));
    }

    @Test
    public void test5TaxonMultiLocus() throws Exception {
        Randomizer.setSeed(1);

        XMLParser parser = new XMLParser();
        beast.core.Runnable runnable = parser.parseFile(
                new File("examples/ACGsimulations/simulateACGs5taxonMultiLocus.xml"));
        runnable.run();

        List<Expectation> expectations = new ArrayList<>();
        expectations.add(new Expectation("acg.CFheight", 1.917, 1e-2));
        expectations.add(new Expectation("acg.CFlength", 4.672, 1e-2));
        expectations.add(new Expectation("acg.nConv", 23.614, 5e-2));

        LogAnalyser logAnalyser = new LogAnalyser("simulateACGs5taxonMultiLocus.stats",
                expectations);

        for (int i=0; i<expectations.size(); i++) {
            assertTrue(expectations.get(i).isValid());
            assertTrue(expectations.get(i).isPassed());
        }

        Files.deleteIfExists(Paths.get("simulateACGs5taxonMultiLocus.stats"));
        Files.deleteIfExists(Paths.get("simulateACGs5taxonMultiLocus.converted"));
        Files.deleteIfExists(Paths.get("simulateACGs5taxonMultiLocus.trees"));
    }
}
