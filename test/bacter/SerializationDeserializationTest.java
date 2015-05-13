/*
 * Copyright (C) 2013 Tim Vaughan <tgvaughan@gmail.com>
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
package bacter;

import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Node;
import beast.util.ClusterTree;
import java.io.ByteArrayInputStream;
import javax.xml.parsers.DocumentBuilderFactory;
import org.junit.Test;
import static org.junit.Assert.*;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;

/**
 * Tests the toString() and fromXML() methods of ConversionGraph.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SerializationDeserializationTest extends TestBase {
    
    public SerializationDeserializationTest() { }
    
    @Test
    public void testXML() throws Exception {
        
        Alignment alignment = getAlignment();
        alignment.setID("alignment");
        
        // ConversionGraph
        ConversionGraph acg = new ConversionGraph();
        ClusterTree tree = new ClusterTree();
        tree.initByName(
                "clusterType", "upgma",
                "taxa", alignment);
        
        acg.assignFrom(tree);
        acg.initByName("alignment", alignment);
        
        //Add recombination event 1
        Node node1 = acg.getExternalNodes().get(0);
        Node node2 = node1.getParent();
        double height1 = 0.5*(node1.getHeight() + node1.getParent().getHeight());
        double height2 = 0.5*(node2.getHeight() + node2.getParent().getHeight());
        int startLocus = 100;
        int endLocus = 200;
        Conversion recomb1 = new Conversion(node1, height1, node2, height2,
                startLocus, endLocus, alignment);
        acg.addConversion(recomb1);
        
        //Add recombination event 2
        node1 = acg.getExternalNodes().get(0);
        node2 = acg.getRoot();
        height1 = 0.5*(node1.getHeight() + node1.getParent().getHeight());
        height2 = node2.getHeight() + 1.0;
        startLocus = 300;
        endLocus = 400;
        Conversion recomb2 = new Conversion(node1, height1, node2, height2,
                startLocus, endLocus, alignment);
        acg.addConversion(recomb2);
        
        // Write ARG out to XML string
        String xmlStr = acg.toXML();
        
        // Build DOM from XML string
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        Document doc = factory.newDocumentBuilder().parse(new ByteArrayInputStream(xmlStr.getBytes()));
        doc.normalize();
        NodeList nodes = doc.getElementsByTagName("*");
        org.w3c.dom.Node docNode = nodes.item(0);
        
        // Read ARG back in from DOM
        ConversionGraph argNew = new ConversionGraph();
        argNew.assignFrom(tree);
        argNew.initByName("alignment", alignment);
        argNew.fromXML(docNode);
        
        // Check that new ARG matches old
        Conversion newRecomb1 = argNew.getConversions(alignment).get(0);
        assertEquals(newRecomb1.getNode1().getNr(),recomb1.getNode1().getNr());
        assertEquals(newRecomb1.getNode2().getNr(),recomb1.getNode2().getNr());
        assertEquals(newRecomb1.getHeight1(),recomb1.getHeight1(), 1e-15);
        assertEquals(newRecomb1.getHeight2(),recomb1.getHeight2(), 1e-15);
        assertEquals(newRecomb1.getStartSite(), recomb1.getStartSite());
        assertEquals(newRecomb1.getEndSite(), recomb1.getEndSite());
        
        Conversion newRecomb2 = argNew.getConversions(alignment).get(1);
        assertEquals(newRecomb2.getNode1().getNr(),recomb2.getNode1().getNr());
        assertEquals(newRecomb2.getNode2().getNr(),recomb2.getNode2().getNr());
        assertEquals(newRecomb2.getHeight1(),recomb2.getHeight1(), 1e-15);
        assertEquals(newRecomb2.getHeight2(),recomb2.getHeight2(), 1e-15);
        assertEquals(newRecomb2.getStartSite(), recomb2.getStartSite());
        assertEquals(newRecomb2.getEndSite(), recomb2.getEndSite());
        
        // Note that there are minor differences in the tree due to
        // rounding errors.  Is this normal!?
    }
    
    @Test
    public void testString() throws Exception {

        Alignment alignment = getAlignment();
        alignment.setID("alignment");
        
        // ConversionGraph
        ConversionGraph arg = new ConversionGraph();
        ClusterTree tree = new ClusterTree();
        tree.initByName(
                "clusterType", "upgma",
                "taxa", alignment);
        
        arg.assignFrom(tree);
        arg.initByName("alignment", alignment);
        
        //Add recombination event 1
        Node node1 = arg.getExternalNodes().get(0);
        Node node2 = node1.getParent();
        double height1 = 0.5*(node1.getHeight() + node1.getParent().getHeight());
        double height2 = 0.5*(node2.getHeight() + node2.getParent().getHeight());
        int startLocus = 100;
        int endLocus = 200;
        Conversion recomb1 = new Conversion(node1, height1, node2, height2,
                startLocus, endLocus, alignment);
        arg.addConversion(recomb1);
        
        //Add recombination event 2
        node1 = arg.getExternalNodes().get(0);
        node2 = arg.getRoot();
        height1 = 0.5*(node1.getHeight() + node1.getParent().getHeight());
        height2 = node2.getHeight() + 1.0;
        startLocus = 300;
        endLocus = 400;
        Conversion recomb2 = new Conversion(node1, height1, node2, height2,
                startLocus, endLocus, alignment);
        arg.addConversion(recomb2);
        
        // Write ARG out to string
        String argString = arg.toString();
        
        // Read ARG back in from string
        ConversionGraph argNew = new ConversionGraph();
        argNew.initByName("alignment", alignment, "fromString", argString);
        
        // Check that new ARG matches old
        Conversion newRecomb1 = argNew.getConversions(alignment).get(0);
        assertEquals(newRecomb1.getNode1().getNr(),recomb1.getNode1().getNr());
        assertEquals(newRecomb1.getNode2().getNr(),recomb1.getNode2().getNr());
        assertEquals(newRecomb1.getHeight1(),recomb1.getHeight1(), 1e-15);
        assertEquals(newRecomb1.getHeight2(),recomb1.getHeight2(), 1e-15);
        assertEquals(newRecomb1.getStartSite(), recomb1.getStartSite());
        assertEquals(newRecomb1.getEndSite(), recomb1.getEndSite());
        
        Conversion newRecomb2 = argNew.getConversions(alignment).get(1);
        assertEquals(newRecomb2.getNode1().getNr(),recomb2.getNode1().getNr());
        assertEquals(newRecomb2.getNode2().getNr(),recomb2.getNode2().getNr());
        assertEquals(newRecomb2.getHeight1(),recomb2.getHeight1(), 1e-15);
        assertEquals(newRecomb2.getHeight2(),recomb2.getHeight2(), 1e-15);
        assertEquals(newRecomb2.getStartSite(), recomb2.getStartSite());
        assertEquals(newRecomb2.getEndSite(), recomb2.getEndSite());
        
        // Note that there are minor differences in the tree due to
        // rounding errors.  Is this normal!?
    }
}