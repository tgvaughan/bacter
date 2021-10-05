package bacter.util;

import bacter.ConversionGraph;
import beast.util.Randomizer;
import beast.util.XMLParser;
import org.junit.Test;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;

import static org.junit.Assert.assertEquals;

public class BacterACGLogReaderTest {

    @Test
    public void test() throws Exception {
        Randomizer.setSeed(1);

        XMLParser parser = new XMLParser();
        beast.core.Runnable runnable = parser.parseFile(
                new File("examples/ACGsimulations/simulateACGs2taxon.xml"));
        runnable.run();

        BacterACGLogReader logReader = new BacterACGLogReader(new File("simulateACGs2taxon.trees"), 0.0);

        assertEquals(100, logReader.getACGCount());

        // The following simply tests whether the ACGs can be loaded without error:
        System.out.println("Root ages:");
        for (ConversionGraph acg : logReader)
            System.out.println(acg.getRoot().getHeight());

        Files.deleteIfExists(Paths.get("simulateACGs2taxon.stats"));
        Files.deleteIfExists(Paths.get("simulateACGs2taxon.converted"));
        Files.deleteIfExists(Paths.get("simulateACGs2taxon.trees"));
    }
}
