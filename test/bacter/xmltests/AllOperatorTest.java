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

import bacter.TestBase;
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

import static org.junit.Assert.assertTrue;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class AllOperatorTest extends TestBase {

    @Test
    public void testSerialSample() throws Exception {
        Randomizer.setSeed(1);

        XMLParser parser = new XMLParser();
        beast.core.Runnable runnable = parser.parseFile(
                new File("examples/allOperatorTests/allOperatorTestSeriallySampled.xml"));
        disableScreenLog(runnable);
        runnable.run();

        List<Expectation> expectations = new ArrayList<>();
        expectations.add(new Expectation("acg.CFheight", 1.909, 0.2));
        expectations.add(new Expectation("acg.CFlength", 4.655, 0.5));
        expectations.add(new Expectation("acg.nConv", 23.381, 0.5));

        LogAnalyser logAnalyser = new LogAnalyser("allOperatorTestSeriallySampled.stats",
                expectations);

        for (Expectation expectation : expectations) {
            assertTrue(expectation.isValid());
            assertTrue(expectation.isPassed());
        }

        Files.deleteIfExists(Paths.get("allOperatorTestSeriallySampled.stats"));
        Files.deleteIfExists(Paths.get("allOperatorTestSeriallySampled.converted"));
        Files.deleteIfExists(Paths.get("allOperatorTestSeriallySampled.trees"));
        Files.deleteIfExists(Paths.get("allOperatorTestSeriallySampled.cf"));
        Files.deleteIfExists(Paths.get("allOperatorTestSeriallySampled.xml.state"));
    }

    @Test
    public void testMultiLocus() throws Exception {
        Randomizer.setSeed(1);

        XMLParser parser = new XMLParser();
        beast.core.Runnable runnable = parser.parseFile(
                new File("examples/allOperatorTests/allOperatorTestSSMultiLocus.xml"));
        disableScreenLog(runnable);
        runnable.run();

        List<Expectation> expectations = new ArrayList<>();
        expectations.add(new Expectation("acg.CFheight", 1.909, 0.2));
        expectations.add(new Expectation("acg.CFlength", 4.655, 0.5));
        expectations.add(new Expectation("acg.nConv", 23.381, 0.5));

        LogAnalyser logAnalyser = new LogAnalyser("allOperatorTestSSMultiLocus.stats",
                expectations);

        for (Expectation expectation : expectations) {
            assertTrue(expectation.isValid());
            assertTrue(expectation.isPassed());
        }

        Files.deleteIfExists(Paths.get("allOperatorTestSSMultiLocus.stats"));
        Files.deleteIfExists(Paths.get("allOperatorTestSSMultiLocus.converted"));
        Files.deleteIfExists(Paths.get("allOperatorTestSSMultiLocus.trees"));
        Files.deleteIfExists(Paths.get("allOperatorTestSSMultiLocus.cf"));
        Files.deleteIfExists(Paths.get("allOperatorTestSSMultiLocus.xml.state"));
    }

    @Test
    public void testShortMultiLocus() throws Exception {
        Randomizer.setSeed(1);

        XMLParser parser = new XMLParser();
        beast.core.Runnable runnable = parser.parseFile(
                new File("examples/allOperatorTests/allOperatorTestShortSSML.xml"));
        disableScreenLog(runnable);
        runnable.run();

        List<Expectation> expectations = new ArrayList<>();
        expectations.add(new Expectation("acg.CFheight", 1.601, 0.2));
        expectations.add(new Expectation("acg.CFlength", 4.17, 0.5));
        expectations.add(new Expectation("acg.nConv", 39.61, 0.5));

        LogAnalyser logAnalyser = new LogAnalyser("allOperatorTestShortSSML.stats",
                expectations);

        for (Expectation expectation : expectations) {
            assertTrue(expectation.isValid());
            assertTrue(expectation.isPassed());
        }

        Files.deleteIfExists(Paths.get("allOperatorTestShortSSML.stats"));
        Files.deleteIfExists(Paths.get("allOperatorTestShortSSML.converted"));
        Files.deleteIfExists(Paths.get("allOperatorTestShortSSML.trees"));
        Files.deleteIfExists(Paths.get("allOperatorTestShortSSML.cf"));
        Files.deleteIfExists(Paths.get("allOperatorTestShortSSML.xml.state"));
    }
}
