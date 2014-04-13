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
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.DataType;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.util.AddOnManager;
import feast.nexus.NexusWriter;
import beast.util.Randomizer;
import beast.util.XMLProducer;
import com.google.common.collect.Lists;
import feast.input.In;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("An alignment produced by simulating sequence evolution down an ARG.")
public class SimulatedAlignment extends Alignment {
    
    public Input<RecombinationGraph> argInput = new In<RecombinationGraph>(
            "arg",
            "Recombination graph down which to simulate evolution.")
            .setRequired();
    
    public Input<SiteModel> siteModelInput = new In<SiteModel>(
            "siteModel",
            "site model for leafs in the beast.tree")
            .setRequired();
    
    public Input<String> outputFileNameInput = In.create("outputFileName",
            "If provided, simulated alignment is additionally written to this file."); 
    
    public Input<Boolean> useNexusInput = new In<Boolean>("useNexus",
            "Use Nexus format to write alignment file.").setDefault(false);
    
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

        // We can't wait for Alignment.initAndValidate() to get the
        // data type for us.
        grabDataType();

        // Simulate alignment
        simulate();
        
        super.initAndValidate();
        
        // Write simulated alignment to disk if requested:
        if (outputFileNameInput.get() != null) {
            PrintStream pstream = new PrintStream(outputFileNameInput.get());
            if (useNexusInput.get())
                NexusWriter.write(this, null, pstream);
            else
                pstream.println(new XMLProducer().toRawXML(this));
        }
    }

    /**
     * Perform actual sequence simulation.
     */
    private void simulate() throws Exception {
        Node cfRoot = arg.getRoot();
        int nTaxa = arg.getLeafNodeCount();
        
        double[] categoryProbs = siteModel.getCategoryProportions(cfRoot);

        int nCategories = siteModel.getCategoryCount();
        int nStates = dataType.getStateCount();
        double[][] transitionProbs = new double[nCategories][nStates*nStates];
        
        int[][] alignment = new int[nTaxa][arg.getSequenceLength()];
        
        for (Recombination recomb : arg.getRecombinations()) {
            int thisLength;
            if (recomb==null)
                thisLength = arg.getClonalFrameSiteCount();
            else
                thisLength = recomb.getSiteCount();
            
            int[] categories = new int[thisLength];
            for (int i=0; i<thisLength; i++)
                categories[i] = Randomizer.randomChoicePDF(categoryProbs);
            
            int[][] regionAlignment = new int[nTaxa][thisLength];
            
            Node thisRoot = arg.getMarginalRoot(recomb);
            
            traverse(recomb, thisRoot, null,
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
            System.out.println(segsites + " segregating sites in region "
                    + "corresponding to recomb " + recomb);
            
            copyToAlignment(alignment, regionAlignment, recomb);
        }
        
        for (int leafIdx=0; leafIdx<nTaxa; leafIdx++) {
            String sSeq = dataType.state2string(alignment[leafIdx]);
            String sTaxon = arg.getNode(leafIdx).getID();
            sequenceInput.setValue(new Sequence(sTaxon, sSeq), this);
        }
    }
    
    /**
     * Traverse a marginal tree simulating a region of the sequence alignment
     * down it.
     * 
     * @param recomb  Recombination identifying the marginal tree
     * @param node Node of the marginal tree
     * @param parentSequence Sequence at the parent node in the marginal tree
     * @param categories Mapping from sites to categories
     * @param transitionProbs
     * @param regionAlignment 
     */
    private void traverse(Recombination recomb, Node node,
            int[] parentSequence,
            int[] categories, double[][] transitionProbs,
            int[][] regionAlignment) {
        
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
            double childHeight = arg.getMarginalNodeHeight(child, recomb);

            // Calculate transition probabilities
            for (int i=0; i<siteModel.getCategoryCount(); i++) {
                siteModel.getSubstitutionModel().getTransitionProbabilities(
                        child, nodeHeight, childHeight,
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
            
            if (arg.isNodeMarginalLeaf(child, recomb)) {
                System.arraycopy(childSequence, 0,
                        regionAlignment[child.getNr()], 0, childSequence.length);
            } else {
                traverse(recomb, child, childSequence,
                        categories, transitionProbs,
                        regionAlignment);
            }
        }
    }
    
    /**
     * Copy region alignment over to full alignment.  Trivial in all
     * cases besides the clonal frame alignment.
     * 
     * @param alignment
     * @param regionAlignment
     * @param recomb 
     */
    private void copyToAlignment(int[][] alignment, int[][] regionAlignment,
            Recombination recomb) {
        
        int nTaxa = alignment.length;
        
        if (recomb==null) {
                
                int ridx=0;
                int sidx=0;
                int jidx=0;
                while (ridx<arg.getNRecombs()) {
                    Recombination nextRecomb = arg.getRecombinations().get(ridx+1);
                    int nextStart = (int)nextRecomb.getStartLocus();
                    if (nextStart>sidx) {
                        for (int leafIdx=0; leafIdx<nTaxa; leafIdx++) {
                            System.arraycopy(regionAlignment[leafIdx], jidx,
                                    alignment[leafIdx], sidx, nextStart-sidx);
                        }
                    }
                    jidx += nextStart-sidx;
                    sidx = (int)nextRecomb.getEndLocus()+1;
                    ridx += 1;
                }
                if (sidx<arg.getSequenceLength()-1) {
                    for (int leafIdx=0; leafIdx<nTaxa; leafIdx++) {
                        System.arraycopy(regionAlignment[leafIdx], jidx,
                                alignment[leafIdx], sidx,
                                arg.getSequenceLength()-sidx);
                    }
                }
                
            } else {
                for (int leafIdx=0; leafIdx<nTaxa; leafIdx++) {
                    System.arraycopy(regionAlignment[leafIdx], 0,
                            alignment[leafIdx], (int)recomb.getStartLocus(),
                            recomb.getSiteCount());
                }
            }
    }
    
    /**
     * HORRIBLE function to identify data type from given description.
     * 
     * @param dataTypeDesc
     * @return DataType instance (null if none found)
     */
    private void grabDataType() {
        if (userDataTypeInput.get() != null) {
            dataType = userDataTypeInput.get();
        } else {

            List<String> dataTypeDescList = Lists.newArrayList();
            List<String> classNames = AddOnManager.find(beast.evolution.datatype.DataType.class, "beast.evolution.datatype");
            for (String className : classNames) {
                try {
                    DataType thisDataType = (DataType) Class.forName(className).newInstance();
                    if (dataTypeInput.get().equals(thisDataType.getDescription())) {
                        dataType = thisDataType;
                        break;
                    }
                    dataTypeDescList.add(thisDataType.getDescription());
                } catch (ClassNotFoundException e) {
                } catch (InstantiationException e) {
                } catch (IllegalAccessException e) {
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
