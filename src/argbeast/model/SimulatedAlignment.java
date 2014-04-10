/*
 * Copyright (C) 2014 Tim Vaughan <tgvaughan@gmail.com>
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

package argbeast.model;

import argbeast.Recombination;
import argbeast.RecombinationGraph;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.Nucleotide;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import feast.input.In;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("An alignment produced by simulating sequence evolution down an ARG.")
public class SimulatedAlignment extends Alignment {
    
    public Input<RecombinationGraph> argInput = new In<RecombinationGraph>(
            "ARG",
            "Recombination graph down which to simulate evolution.")
            .setRequired();
    
    public Input<SiteModel> siteModelInput = new In<SiteModel>(
            "siteModel",
            "site model for leafs in the beast.tree")
            .setRequired();
    
    public Input<DataType> dataTypeInput = new In<DataType>(
            "dataType", "Data type for alignment").setDefault(new Nucleotide());
    
    public Input<String> outputFileNameInput = In.create("outputFileName",
            "If provided, simulated alignment is additionally written to this file.");    
    
    private RecombinationGraph arg;
    private SiteModel siteModel;
    private DataType dataType;
    
    public SimulatedAlignment() {
        In.setOptional(sequenceInput);
    }

    @Override
    public void initAndValidate() throws Exception {

        arg = argInput.get();
        siteModel = siteModelInput.get();
        dataType = dataTypeInput.get();
        
        simulate();
        
        super.initAndValidate();
    }

    /**
     * Perform actual sequence simulation.
     */
    private void simulate() {
        RecombinationGraph arg = argInput.get();
        Node cfRoot = arg.getRoot();
        int nTaxa = arg.getLeafNodeCount();
        
        double[] categoryProbs = siteModel.getCategoryProportions(cfRoot);
        int[] categories = new int[arg.getSequenceLength()];
        for (int i = 0; i < arg.getSequenceLength(); i++) {
            categories[i] = Randomizer.randomChoicePDF(categoryProbs);
        }
        int nCategories = siteModel.getCategoryCount();
        int nStates = dataType.getStateCount();
        double[][] transitionProbs = new double[nCategories][nStates*nStates];
        
        int[][] alignment = new int[nTaxa][arg.getSequenceLength()];
        
        for (Recombination recomb : arg.getRecombinations()) {

        }
    }
    
    
    private void traverse(Node node, int[] parentSequence,
            int[] categories, double[][] transitionProbs,
            Recombination recomb) {
        
        // Draw random ancestral sequence if necessary
        if (parentSequence == null) {
            if (recomb == null)
                parentSequence = new int[arg.getClonalFrameSiteCount()];
            else
                parentSequence = new int[recomb.getSiteCount()];
            
            double[] frequencies = siteModel.getSubstitutionModel().getFrequencies();
            for (int i=0; i<parentSequence.length; i++)
                parentSequence[i] = Randomizer.randomChoicePDF(frequencies);
        }
        
        double nodeHeight = arg.getMarginalNodeHeight(node, recomb);
        for (Node child : arg.getMarginalChildren(node, recomb)) {
            double childHeight = arg.getMarginalNodeHeight(node, recomb);

            // Calculate transition probabilities
            for (int i=0; i<siteModel.getCategoryCount(); i++) {
                siteModel.getSubstitutionModel().getTransitionProbabilities(
                        child, nodeHeight, childHeight,
                        siteModel.getRateForCategory(i, child),
                        transitionProbs[i]);
            }
        }
    }
}
