package bacter.model;

import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.Locus;
import bacter.TestBase;
import bacter.acgannotator.ACGAnnotator;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.ClusterTree;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

public class ACGLikelihoodCircularVsLinearTest extends TestBase {

    public ACGLikelihoodCircularVsLinearTest() { }

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
        int startLocus = 6; //9;
        int endLocus = 11; //2;

        Conversion recomb1 = new Conversion(node1, height1, node2, height2,
                startLocus, endLocus, acg, locus);
        acg.addConversion(recomb1);

        double logP = argLikelihood.calculateLogP();


        ConversionGraph acg2 = new ConversionGraph();
        acg2.assignFrom(tree);
        acg2.initByName("locus", locus,
                "circularGenome", "false",
                "endSiteBetaBinom", "false");

        // Likelihood

        ACGLikelihood argLikelihood2 = new ACGLikelihood();
        argLikelihood2.initByName(
                "locus", locus,
                "tree", acg2,
                "siteModel", siteModel);

        acg2.setEverythingDirty(true);

        //Add a single recombination event
        acg2.getNodesAsArray();
        Node node21 = acg2.getExternalNodes().get(0); //acg.getExternalNodes().get(0);
        Node node22 = node21.getParent();
        double height21 = 0.5*(node21.getHeight() + node21.getParent().getHeight());
        double height22 = 0.5*(node22.getHeight() + node22.getParent().getHeight());
        int startLocus2 = 6; //9;
        int endLocus2 = 11; //2;

        Conversion recomb21 = new Conversion(node21, height21, node22, height22,
                startLocus2, endLocus2, acg2, locus);
        acg2.addConversion(recomb21);

        double logP2 = argLikelihood2.calculateLogP();

        double relativeDiff2 = Math.abs(2.0*(logP2-logP)/(logP2+logP));

        assertTrue(relativeDiff2<1e-14);

    }

    }
