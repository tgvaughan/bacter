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

import bacter.ConversionGraph;
import bacter.Locus;
import bacter.TestBase;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.util.ClusterTree;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ACGLikelihoodApproxTest extends TestBase {

    @Test
    public void testPairwiseDistances() throws Exception {

        List<Sequence> sequenceList = new ArrayList<>();
                                           //01234567890123456789
        sequenceList.add(new Sequence("t1", "AAAAAAAAAAAAAAAAAAAA"));
        sequenceList.add(new Sequence("t2", "AAAACAAAAAAGAAAAAAAA"));
        sequenceList.add(new Sequence("t3", "CCCCCCCCCCCCCCCCCCCC"));
        Alignment alignment = new Alignment(sequenceList, "nucleotide");

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

        Assert.assertEquals(likelihoodApprox.getPairwiseDistance(0, 1, 0, 4), 0);
        Assert.assertEquals(likelihoodApprox.getPairwiseDistance(0,1,0,5), 1);
        Assert.assertEquals(likelihoodApprox.getPairwiseDistance(0,1,0,11), 1);
        Assert.assertEquals(likelihoodApprox.getPairwiseDistance(0,1,0,12), 2);
        Assert.assertEquals(likelihoodApprox.getPairwiseDistance(0,1,0,20), 2);

        Assert.assertEquals(likelihoodApprox.getPairwiseDistance(1,2,2,8), 5);
        Assert.assertEquals(likelihoodApprox.getPairwiseDistance(1,2,0,20), 19);
        Assert.assertEquals(likelihoodApprox.getPairwiseDistance(0,2,2,8), 6);
        Assert.assertEquals(likelihoodApprox.getPairwiseDistance(0,2,0,20), 20);
    }
}
