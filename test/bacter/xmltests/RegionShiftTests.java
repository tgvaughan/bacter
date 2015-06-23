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
public class RegionShiftTests {

    @Test
    public void testRegionShift() throws Exception {
        Randomizer.setSeed(1);

        XMLParser parser = new XMLParser();
        beast.core.Runnable runnable = parser.parseFile(
                new File("examples/regionShiftTests/regionShiftTest5taxon.xml"));
        runnable.run();

        List<Expectation> expectations = new ArrayList<>();
        expectations.add(new Expectation("acg.CFheight", 1.606, 0.2));
        expectations.add(new Expectation("acg.CFlength", 4.181, 0.5));
        expectations.add(new Expectation("acg.nConv", 21.0, 0.5));

        LogAnalyser logAnalyser = new LogAnalyser("regionShiftTest5taxon.stats",
                expectations);

        for (int i=0; i<expectations.size(); i++) {
            assertTrue(expectations.get(i).isValid());
            assertTrue(expectations.get(i).isPassed());
        }

        Files.deleteIfExists(Paths.get("regionShiftTest5taxon.stats"));
        Files.deleteIfExists(Paths.get("regionShiftTest5taxon.converted"));
        Files.deleteIfExists(Paths.get("regionShiftTest5taxon.trees"));
        Files.deleteIfExists(Paths.get("regionShiftTest5taxon.cf"));
        Files.deleteIfExists(Paths.get("regionShiftTest5taxon.xml.state"));
    }

    @Test
    public void testBoundaryShift() throws Exception {
        Randomizer.setSeed(1);

        XMLParser parser = new XMLParser();
        beast.core.Runnable runnable = parser.parseFile(
                new File("examples/regionShiftTests/boundaryShiftTest5taxon.xml"));
        runnable.run();

        List<Expectation> expectations = new ArrayList<>();
        expectations.add(new Expectation("acg.CFheight", 1.606, 0.2));
        expectations.add(new Expectation("acg.CFlength", 4.181, 0.5));
        expectations.add(new Expectation("acg.nConv", 21.0, 0.5));

        LogAnalyser logAnalyser = new LogAnalyser("boundaryShiftTest5taxon.stats",
                expectations);

        for (int i=0; i<expectations.size(); i++) {
            assertTrue(expectations.get(i).isValid());
            assertTrue(expectations.get(i).isPassed());
        }

        Files.deleteIfExists(Paths.get("boundaryShiftTest5taxon.stats"));
        Files.deleteIfExists(Paths.get("boundaryShiftTest5taxon.converted"));
        Files.deleteIfExists(Paths.get("boundaryShiftTest5taxon.trees"));
        Files.deleteIfExists(Paths.get("boundaryShiftTest5taxon.cf"));
        Files.deleteIfExists(Paths.get("boundaryShiftTest5taxon.xml.state"));
    }

}
