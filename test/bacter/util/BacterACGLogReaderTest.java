package bacter.util;

import bacter.ConversionGraph;
import beast.util.Randomizer;
import beast.util.XMLParser;
import beast.util.XMLParserException;
import org.junit.Test;
import org.xml.sax.SAXException;
import test.beast.beast2vs1.trace.Expectation;
import test.beast.beast2vs1.trace.LogAnalyser;

import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class BacterACGLogReaderTest {

    @Test
    public void test() throws Exception {

        Randomizer.setSeed(1);

        XMLParser parser = new XMLParser();
        beast.core.Runnable runnable = parser.parseFile(
                new File("examples/ACGsimulations/simulateACGs2taxon.xml"));
        runnable.run();

        BacterACGLogReader logReader = new BacterACGLogReader(new File("simulateACGs2taxon.trees"), 0.0);

        System.out.println("ACG count: " + logReader.getACGCount());

        System.out.println("Root ages:");
        for (ConversionGraph acg : logReader)
            System.out.println(acg.getRoot().getHeight());

        Files.deleteIfExists(Paths.get("simulateACGs2taxon.stats"));
        Files.deleteIfExists(Paths.get("simulateACGs2taxon.converted"));
        Files.deleteIfExists(Paths.get("simulateACGs2taxon.trees"));
    }
}
