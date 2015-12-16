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

package bacter.model;

import bacter.*;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.util.ClusterTree;
import beast.util.TreeParser;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ACGLikelihoodApproxTest extends TestBase {

    @Test
    public void testPairwiseDistances() throws Exception {

        List<Sequence> sequences = new ArrayList<>();
                                        //01234567890123456789
        sequences.add(new Sequence("t1", "AAAAAAAAAAAAAAAAAAAA"));
        sequences.add(new Sequence("t2", "AAAACAAAAAAGAAAAAAAA"));
        sequences.add(new Sequence("t3", "CCCCCCCCCCCCCCCCCCCC"));
        Alignment alignment = new Alignment(sequences, "nucleotide");

        Locus locus = new Locus("locus", alignment);
        ConversionGraph acg = new ConversionGraph();
        ClusterTree tree = new ClusterTree();
        tree.initByName(
                "clusterType", "upgma",
                "taxa", locus.getAlignment());

        acg.assignFrom(tree);
        acg.initByName("locus", locus);

        ACGLikelihoodApprox likelihoodApprox = new ACGLikelihoodApprox();
        likelihoodApprox.initByName(
                "acg", acg,
                "substitutionRate", "1.0",
                "alignment", alignment,
                "locus", locus);

        Assert.assertEquals(0, likelihoodApprox.getPairwiseDistance(0, 1, 0, 4));
        Assert.assertEquals(1, likelihoodApprox.getPairwiseDistance(0,1,0,5));
        Assert.assertEquals(1, likelihoodApprox.getPairwiseDistance(0,1,0,11));
        Assert.assertEquals(2, likelihoodApprox.getPairwiseDistance(0,1,0,12));
        Assert.assertEquals(2, likelihoodApprox.getPairwiseDistance(0,1,0,20));

        Assert.assertEquals(5, likelihoodApprox.getPairwiseDistance(1,2,2,8));
        Assert.assertEquals(19, likelihoodApprox.getPairwiseDistance(1,2,0,20));
        Assert.assertEquals(6, likelihoodApprox.getPairwiseDistance(0,2,2,8));
        Assert.assertEquals(20, likelihoodApprox.getPairwiseDistance(0,2,0,20));
    }

    @Test
    public void testTreeHeightMap() throws Exception {

        List<Sequence> sequences = new ArrayList<>();
                                        //01234567890123456789
        sequences.add(new Sequence("t1", "GGGGGGGGGGGGGGGGGGGG"));
        sequences.add(new Sequence("t2", "CCCCCCCCCCCCCCCCCCCC"));
        sequences.add(new Sequence("t3", "TTTTTTTTTTTTTTTTTTTT"));
        Alignment alignment = new Alignment(sequences, "nucleotide");
        Locus locus = new Locus("locus", alignment);

        TreeParser tree = new TreeParser(alignment, "((t1:1,t2:1):1,t3:2):0;");
        ConversionGraph acg = new ConversionGraph();
        acg.assignFrom(tree);
        acg.initByName("locus", locus);

        ACGLikelihoodApprox likelihoodApprox = new ACGLikelihoodApprox();
        likelihoodApprox.initByName(
                "acg", acg,
                "substitutionRate", "1.0",
                "alignment", alignment,
                "locus", locus);

        Map<Double, Coalescence> heightMap = likelihoodApprox.getCoalescenceHeights();

        Assert.assertEquals(2, heightMap.size());
        Assert.assertTrue(heightMap.containsKey(1.0));
        Assert.assertTrue(heightMap.containsKey(2.0));
        Assert.assertTrue(heightMap.get(1.0).equals(new Coalescence("[0,20]{0}{1}")));
        Assert.assertTrue(heightMap.get(2.0).equals(new Coalescence("[0,20]{0,1}{2}")));
    }

    @Test
    public void testACGHeightMap() throws Exception {

        List<Sequence> sequences = new ArrayList<>();
        //01234567890123456789
        sequences.add(new Sequence("t1", "GGGGGGGGGGGGGGGGGGGG"));
        sequences.add(new Sequence("t2", "CCCCCCCCCCCCCCCCCCCC"));
        sequences.add(new Sequence("t3", "TTTTTTTTTTTTTTTTTTTT"));
        Alignment alignment = new Alignment(sequences, "nucleotide");
        Locus locus = new Locus("locus", alignment);

        TreeParser tree = new TreeParser(alignment, "((t1:1,t2:1):1,t3:2):0;");
        ConversionGraph acg = new ConversionGraph();
        acg.assignFrom(tree);
        acg.initByName("locus", locus);

        Conversion conversion = new Conversion();
        conversion.setNode1(acg.getNode(0));
        conversion.setHeight1(0.5);
        conversion.setNode2(acg.getNode(2));
        conversion.setHeight2(1.5);
        conversion.setStartSite(0);
        conversion.setEndSite(9);
        conversion.setLocus(locus);
        acg.addConversion(conversion);

        System.out.println(acg);

        ACGLikelihoodApprox likelihoodApprox = new ACGLikelihoodApprox();
        likelihoodApprox.initByName(
                "acg", acg,
                "substitutionRate", "1.0",
                "alignment", alignment,
                "locus", locus);

        Map<Double, Coalescence> heightMap = likelihoodApprox.getCoalescenceHeights();

        Assert.assertEquals(3, heightMap.size());
        Assert.assertTrue(heightMap.containsKey(1.0));
        Assert.assertTrue(heightMap.containsKey(1.5));
        Assert.assertTrue(heightMap.containsKey(2.0));
        Assert.assertTrue("height map contains incorrect coalescence.",
                heightMap.get(1.0).equals(new Coalescence("[10,20]{0}{1}")));
        Assert.assertTrue("height map contains incorrect coalescence.",
                heightMap.get(1.5).equals(new Coalescence("[0,10]{0}{2}")));
        Assert.assertTrue("height map contains incorrect coalescence.",
                heightMap.get(2.0).equals(new Coalescence("[0,10]{0,2}{1} [10,20]{0,1}{2}")));
    }
}
