package bacter.model;

import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.Locus;
import bacter.TestBase;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.evolution.tree.ClusterTree;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

public class ACGCoalescentCircularTest  extends TestBase {

    public ACGCoalescentCircularTest() { }

    @Test
    public void testLikelihoodFixedData() throws Exception {

        Locus locus = new Locus("locus", getSimpleAlignment());

        // ConversionGraph
        ConversionGraph acg = new ConversionGraph();
        ClusterTree tree = new ClusterTree();
        tree.initByName(
                "clusterType", "upgma",
                "taxa", locus.getAlignment());

        acg.assignFrom(tree);
        acg.initByName("locus", locus,
                                "circularGenome", "true",
                                "endSiteBetaBinom", "false");


        //Add a single recombination event
        acg.getNodesAsArray();
        Node node1 = acg.getExternalNodes().get(0); //acg.getExternalNodes().get(0);
        Node node2 = node1.getParent();
        double height1 = 0.5*(node1.getHeight() + node1.getParent().getHeight());
        double height2 = 0.5*(node2.getHeight() + node2.getParent().getHeight());
        int startLocus = 9;
        int endLocus = 2;
        Conversion recomb1 = new Conversion(node1, height1, node2, height2,
                startLocus, endLocus, acg, locus);
        acg.addConversion(recomb1);


        ConstantPopulation popFunc = new ConstantPopulation();
        popFunc.initByName("popSize", new RealParameter("1.0"));

        ACGCoalescent acgCoalescent = new ACGCoalescent();

        acgCoalescent.initByName(
                "populationModel", popFunc, "rho", "1e-6", "delta", "3", "tree", acg, "circularGenome", "true");

        double logP = acgCoalescent.calculateLogP();


        //Add a single recombination event
        startLocus = 3;
        endLocus = 8;
        Conversion recomb2 = new Conversion(node1, height1, node2, height2,
                startLocus, endLocus, acg, locus);
        acg.deleteConversion(recomb1);
        acg.addConversion(recomb2);

        ACGCoalescent acgCoalescent2 = new ACGCoalescent();

        RealParameter delta = new RealParameter();
        delta.initByName("value", "3.", "upper","100");

        acgCoalescent2.initByName(
                "populationModel", popFunc, "rho", "1e-6", "delta", delta, "tree", acg,"circularGenome", "true");

        double logPOtherHalf = acgCoalescent2.calculateLogP();

        double relativeDiff = Math.abs(2.0*(logPOtherHalf-logP)/(logPOtherHalf+logP));

        assertTrue(relativeDiff<1e-14);
    }

    }
