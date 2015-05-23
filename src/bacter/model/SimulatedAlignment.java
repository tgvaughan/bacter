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

package bacter.model;

import bacter.ConversionGraph;
import bacter.Locus;
import bacter.MarginalTree;
import bacter.Region;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.DataType;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.util.AddOnManager;
import beast.util.Randomizer;
import beast.util.XMLProducer;
import feast.nexus.CharactersBlock;
import feast.nexus.NexusBuilder;
import feast.nexus.TaxaBlock;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("An alignment produced by simulating sequence evolution down an ACG.")
public class SimulatedAlignment extends Alignment {
    
    public Input<ConversionGraph> acgInput = new Input<>(
            "acg",
            "Conversions graph down which to simulate evolution.",
            Input.Validate.REQUIRED);
    
    public Input<SiteModel> siteModelInput = new Input<>(
            "siteModel",
            "site model for leafs in the beast.tree",
            Input.Validate.REQUIRED);
    
    public Input<String> outputFileNameInput = new Input<>(
            "outputFileName",
            "If provided, simulated alignment is additionally written to this file."); 
    
    public Input<Boolean> useNexusInput = new Input<>(
            "useNexus",
            "Use Nexus format to write alignment file.",
            false);

    public Input<Locus> locusInput = new Input<>("locus",
            "Locus for which alignment will be simulated.");
    
    private ConversionGraph acg;
    private SiteModel siteModel;
    private DataType dataType;
    
    public SimulatedAlignment() {
        sequenceInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() throws Exception {

        acg = acgInput.get();
        siteModel = siteModelInput.get();

        Locus locus;
        if (acg.getLoci().size()==1)
            locus = acg.getLoci().get(0);
        else {
            locus = locusInput.get();
            if (locus == null)
                throw new IllegalArgumentException("Must specify locus" +
                        " when simulating alignment from ACG associated" +
                        " with multiple loci.");
        }

        // We can't wait for Alignment.initAndValidate() to get the
        // data type for us.
        grabDataType();

        // Simulate alignment
        simulate(locus);
        
        super.initAndValidate();
        
        // Write simulated alignment to disk if requested:
        if (outputFileNameInput.get() != null) {
            PrintStream pstream = new PrintStream(outputFileNameInput.get());
            if (useNexusInput.get()) {
                NexusBuilder nb = new NexusBuilder();
                nb.append(new TaxaBlock(acg.getTaxonset()));
                nb.append(new CharactersBlock(this));
                nb.write(pstream);
            } else
                pstream.println(new XMLProducer().toRawXML(this));
        }
    }

    /**
     * Perform actual sequence simulation.
     */
    private void simulate(Locus locus) throws Exception {
        Node cfRoot = acg.getRoot();
        int nTaxa = acg.getLeafNodeCount();
        
        double[] categoryProbs = siteModel.getCategoryProportions(cfRoot);

        int nCategories = siteModel.getCategoryCount();
        int nStates = dataType.getStateCount();
        double[][] transitionProbs = new double[nCategories][nStates*nStates];
        
        int[][] alignment = new int[nTaxa][locus.getSiteCount()];
        
        for (Region region : acg.getRegions(locus)) {
            int thisLength = region.getRegionLength();

            MarginalTree marginalTree = new MarginalTree(acg, region);
            
            int[] categories = new int[thisLength];
            for (int i=0; i<thisLength; i++)
                categories[i] = Randomizer.randomChoicePDF(categoryProbs);
            
            int[][] regionAlignment = new int[nTaxa][thisLength];
            
            Node thisRoot = marginalTree.getRoot();
            
            int[] parentSequence = new int[region.getRegionLength()];
            double[] frequencies = siteModel.getSubstitutionModel().getFrequencies();
            for (int i=0; i<parentSequence.length; i++)
                parentSequence[i] = Randomizer.randomChoicePDF(frequencies);

            traverse(thisRoot, parentSequence,
                    categories, transitionProbs,
                    regionAlignment);
            
            // DEBUG: Count differences
            int segsites = 0;
            for (int s=0; s<regionAlignment[0].length; s++) {
                int state = regionAlignment[0][s];
                for (int l=1; l<nTaxa; l++) {
                    if (state != regionAlignment[l][s]) {
                        segsites += 1;
                        break;
                    }
                }
            }
            System.out.println(segsites + " segregating sites in region " + region);
            
            copyToAlignment(alignment, regionAlignment, region);
        }
        
        for (int leafIdx=0; leafIdx<nTaxa; leafIdx++) {
            String sSeq = dataType.state2string(alignment[leafIdx]);
            String sTaxon = acg.getNode(leafIdx).getID();
            sequenceInput.setValue(new Sequence(sTaxon, sSeq), this);
        }
    }
    
    /**
     * Traverse a marginal tree simulating a region of the sequence alignment
     * down it.
     * 
     * @param node Node of the marginal tree
     * @param parentSequence Sequence at the parent node in the marginal tree
     * @param categories Mapping from sites to categories
     * @param transitionProbs
     * @param regionAlignment 
     */
    private void traverse(Node node,
            int[] parentSequence,
            int[] categories, double[][] transitionProbs,
            int[][] regionAlignment) {
        
        for (Node child : node.getChildren()) {

            // Calculate transition probabilities
            for (int i=0; i<siteModel.getCategoryCount(); i++) {
                siteModel.getSubstitutionModel().getTransitionProbabilities(
                        child, node.getHeight(), child.getHeight(),
                        siteModel.getRateForCategory(i, child),
                        transitionProbs[i]);
            }
            
            // Draw characters on child sequence
            int[] childSequence = new int[parentSequence.length];
            int nStates = dataType.getStateCount();
            double[] charProb = new double[nStates];
            for (int i=0; i<childSequence.length; i++) {
                int category = categories[i];
                System.arraycopy(transitionProbs[category],
                        parentSequence[i]*nStates, charProb, 0, nStates);
                childSequence[i] = Randomizer.randomChoicePDF(charProb);
            }
            
            if (child.isLeaf()) {
                System.arraycopy(childSequence, 0,
                        regionAlignment[child.getNr()], 0, childSequence.length);
            } else {
                traverse(child, childSequence,
                        categories, transitionProbs,
                        regionAlignment);
            }
        }
    }
    
    /**
     * Copy region alignment over to full alignment.
     *
     * @param alignment
     * @param regionAlignment
     * @param region
     */
    private void copyToAlignment(int[][] alignment, int[][] regionAlignment,
            Region region) {
        
        int nTaxa = alignment.length;
        
        for (int leafIdx=0; leafIdx<nTaxa; leafIdx++) {
            System.arraycopy(regionAlignment[leafIdx], 0,
                    alignment[leafIdx], region.leftBoundary,
                    region.getRegionLength());
        }
    }
    
    /**
     * HORRIBLE function to identify data type from given description.
     * 
     * @return DataType instance (null if none found)
     */
    private void grabDataType() {
        if (userDataTypeInput.get() != null) {
            dataType = userDataTypeInput.get();
        } else {

            List<String> dataTypeDescList = new ArrayList<>();
            List<String> classNames = AddOnManager.find(beast.evolution.datatype.DataType.class, "beast.evolution.datatype");
            for (String className : classNames) {
                try {
                    DataType thisDataType = (DataType) Class.forName(className).newInstance();
                    if (dataTypeInput.get().equals(thisDataType.getTypeDescription())) {
                        dataType = thisDataType;
                        break;
                    }
                    dataTypeDescList.add(thisDataType.getTypeDescription());
                } catch (ClassNotFoundException
                    | InstantiationException
                    | IllegalAccessException e) {
                }
            }
            if (dataType == null) {
                throw new IllegalArgumentException("Data type + '"
                        + dataTypeInput.get()
                        + "' cannot be found.  Choose one of "
                        + Arrays.toString(dataTypeDescList.toArray(new String[0])));
            }
        }
    }
}
