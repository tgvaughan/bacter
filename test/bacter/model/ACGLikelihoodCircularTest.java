package bacter.model;

import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.Locus;
import bacter.TestBase;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.tree.Node;
//import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.base.evolution.tree.ClusterTree;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

public class ACGLikelihoodCircularTest extends TestBase {

    public ACGLikelihoodCircularTest() { }

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
                                "endSiteBetaBinom", "true");

        // Site model:
        JukesCantor jc = new JukesCantor();
        jc.initByName();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName(
                "substModel", jc);

        // Likelihood

        ACGLikelihood argLikelihood = new ACGLikelihood();
        argLikelihood.initByName(
                "locus", locus,
                "tree", acg,
                "siteModel", siteModel);

        acg.setEverythingDirty(true);

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

        //new
        Conversion recomb12 = new Conversion(node1, height1, node2, height2,
                0, 1, acg, locus);
        acg.addConversion(recomb12);

        double logP = argLikelihood.calculateLogP();

        //Define recombination on other half of circular genome
        acg.deleteConversion(recomb1);
        acg.deleteConversion(recomb12);

        startLocus = 3; //7;
        endLocus = 8; //13;
        Conversion recomb2 = new Conversion(node1, height1, node2, height2,
                startLocus, endLocus, acg, locus);
        acg.addConversion(recomb2);
        //acg.initAndValidate();
        //argLikelihood.initAndValidate();

        //new
        Conversion recomb22 = new Conversion(node1, height1, node2, height2,
                6, 7, acg, locus);
        acg.addConversion(recomb22);

        double logPOtherHalf = argLikelihood.calculateLogP();

        double relativeDiff = Math.abs(2.0*(logPOtherHalf-logP)/(logPOtherHalf+logP));

        assertTrue(relativeDiff<1e-14);

    }

    }
