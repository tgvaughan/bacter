package bacter.model;

import bacter.ConversionGraph;
import beast.evolution.likelihood.GenericTreeLikelihood;

import java.util.HashSet;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ACGLikelihoodNew extends GenericTreeLikelihood {

    HashSet<LocalNode> marginalNodes;
    List<LocalNode> marginalRoots;

    ConversionGraph acg;

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        if (!(treeInput.get() instanceof ConversionGraph))
            throw new IllegalArgumentException("Tree input of ACGLikelihood " +
                    "must be an instance of ConversionGraph.");

        acg = (ConversionGraph)treeInput.get();
    }

    public void buildMarginalTrees() {


    }

}
