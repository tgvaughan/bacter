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
import beagle.Beagle;
import beagle.BeagleFactory;
import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.BeerLikelihoodCore;
import beast.evolution.likelihood.BeerLikelihoodCore4;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.likelihood.LikelihoodCore;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import com.google.common.collect.LinkedHashMultiset;
import com.google.common.collect.Multiset;

import java.util.*;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Probability of sequence data given recombination graph.")
public class ACGLikelihoodBeagle extends GenericTreeLikelihood {

    public Input<Locus> locusInput = new Input<>(
            "locus",
            "Locus associated with alignment to evaluate probability of.",
            Input.Validate.REQUIRED);

    public Input<Boolean> useAmbiguitiesInput = new Input<>(
            "useAmbiguities",
            "Whether sites containing ambiguous states should be handled " +
                    "instead of ignored (the default)", false);

    protected ConversionGraph acg;

    protected SiteModel.Base siteModel;
    protected BranchRateModel branchRateModel;
    protected SubstitutionModel.Base substitutionModel;
    protected Alignment alignment;
    protected Locus locus;
    protected int nStates;

    protected Map<Region, Multiset<int[]>> patterns;
    protected Map<Region, Multiset<int[]>> storedPatterns;
    protected Map<Region, List<Integer>> constantPatterns;
    protected Map<Region, List<Integer>> storedConstantPatterns;
    protected Map<Region, Beagle> beagleInstances;
    protected Map<Region, Double> regionLogLikelihoods;
    protected Map<Region, Double> storedRegionLogLikelihoods;

    int[] nodeNrs;
    double[] edgeLengths;
    int[] operationList, operationListIdx;

    public ACGLikelihoodBeagle() {
        // We allow alignments to be specified using Locus objects.
        dataInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() throws Exception {

        if (treeInput.get() instanceof ConversionGraph)
            acg = (ConversionGraph)treeInput.get();
        else
            throw new IllegalArgumentException("'Tree' input to ACGLikelihood must " +
                    "be of type ConversionGraph.");

        locus = locusInput.get();
        if (locus.hasAlignment()) {
            alignment = locus.getAlignment();
        } else {
            if (dataInput.get() != null)
                alignment = dataInput.get();
            else
                throw new IllegalArgumentException("No alignment associated with " +
                        "locus " + locus.getID() + " provided to ACGLikelihood " +
                        "and none given explicitly.");
        }

        nStates = alignment.getMaxStateCount();

        siteModel = (SiteModel.Base) siteModelInput.get();
        substitutionModel = (SubstitutionModel.Base) siteModel.getSubstitutionModel();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();

            if (!(branchRateModel instanceof StrictClockModel))
                throw new IllegalArgumentException("ACGLikelihood currently only" +
                        "supports strict clock models.");
        } else
            branchRateModel = new StrictClockModel();

        patterns = new HashMap<>();
        storedPatterns = new HashMap<>();
        constantPatterns = new HashMap<>();
        storedConstantPatterns = new HashMap<>();
        beagleInstances = new HashMap<>();
        regionLogLikelihoods = new HashMap<>();
        storedRegionLogLikelihoods = new HashMap<>();

        edgeLengths = new double[acg.getNodeCount()-1];
        nodeNrs = new int[acg.getNodeCount()-1];
        for (int i=0; i<nodeNrs.length; i++)
            nodeNrs[i] = i;

        operationList = new int[acg.getInternalNodeCount() * Beagle.OPERATION_TUPLE_SIZE];
        operationListIdx = new int[1];
    }

    @Override
    public double calculateLogP() {
        updatePatterns();
        updateBeagleInstances();

        logP = 0.0;

        regionLogLikelihoods.keySet().retainAll(acg.getRegions(locus));

        int rootNr = acg.getRoot().getNr();
        double[] regionLogP = new double[1];

        for (Region region : acg.getRegions(locus)) {

            if (!regionLogLikelihoods.containsKey(region)) {
                Beagle beagle = beagleInstances.get(region);
                MarginalTree marginalTree = new MarginalTree(acg, region);
                traverse(beagle, marginalTree.getRoot(), region);

                beagle.setCategoryRates(siteModel.getCategoryRates(null));
                beagle.setCategoryWeights(0, siteModel.getCategoryProportions(null));
                beagle.setStateFrequencies(0, substitutionModel.getFrequencies());

                beagle.updateTransitionMatrices(0, nodeNrs,
                        null, null, edgeLengths, edgeLengths.length);

                beagle.updatePartials(operationList, operationList.length, Beagle.NONE);

                beagle.calculateRootLogLikelihoods(
                        new int[]{rootNr},
                        new int[]{0},
                        new int[]{0},
                        new int[]{Beagle.NONE},
                        1,
                        regionLogP);

                regionLogLikelihoods.put(region, regionLogP[0]);
                logP += regionLogP[0];
            } else {
                logP += regionLogLikelihoods.get(region);
            }
        }

        return logP;
    }

    /**
     * Ensure pattern counts are up to date.
     */
    private void updatePatterns() {
        List<Region> regionList = acg.getRegions(locus);

        // Remove stale pattern sets
        patterns.keySet().retainAll(regionList);
        constantPatterns.keySet().retainAll(regionList);

        for (Region region : regionList) {

            if (patterns.containsKey(region))
                continue;

            // Add new pattern set
            Multiset<int[]> patSet = LinkedHashMultiset.create();
            for (int j=region.leftBoundary; j<region.rightBoundary; j++) {
                int [] pat = alignment.getPattern(alignment.getPatternIndex(j));
                patSet.add(pat);
            }
            patterns.put(region, patSet);

            // Compute corresponding constant pattern list
            List<Integer> constantPatternList = new ArrayList<>();

            int patternIdx = 0;
            for (int[] pattern : patSet.elementSet()) {
                boolean isConstant = true;
                for (int i=1; i<pattern.length; i++)
                    if (pattern[i] != pattern[0]) {
                        isConstant = false;
                        break;
                    }

                if (isConstant) {
                    if (alignment.getDataType().isAmbiguousState(pattern[0])) {
                        if (useAmbiguitiesInput.get()) {
                            for (int state : alignment.getDataType().getStatesForCode(pattern[0]))
                                constantPatternList.add(patternIdx * nStates + state);
                        }
                    } else {
                        constantPatternList.add(patternIdx * nStates + pattern[0]);
                    }
                }

                patternIdx += 1;
            }

            constantPatterns.put(region, constantPatternList);
        }
    }
    
    
    /**
     * Initialize beagle instances.
     */
    private void updateBeagleInstances() {

        List<Region> regionList = acg.getRegions(locus);
        beagleInstances.keySet().retainAll(regionList);

        for (Region region : regionList) {
            Beagle beagleInstance = BeagleFactory.loadBeagleInstance(
                    acg.getLeafNodeCount(), // Number of tips
                    acg.getNodeCount(), // Number of partials
                    useAmbiguitiesInput.get() ? 0 : acg.getLeafNodeCount(), // Number of compacts
                    nStates, // Number of discrete states in model (4 for DNA)
                    patterns.get(region).elementSet().size(), // Number of patterns
                    1, // Number of eigen decompositions
                    acg.getNodeCount()-1, // Number of transition matrices (one per edge)
                    siteModel.getCategoryCount(), // Number of rate categories
                    0, // Number of scaling buffers (0 means not needed)
                    null, // Potential resource list (null -> no restriction)
                    0, // bit flags indicating preferred implementation characteristics
                    0); // bit flags indicating required implementation characteristics

            if (useAmbiguitiesInput.get()) {
                setPartials(beagleInstance, patterns.get(region));
            } else {
                setStates(beagleInstance, patterns.get(region));
            }

            double weights[] = new double[patterns.get(region).elementSet().size()];
            int i=0;
            for (int[] pattern : patterns.get(region).elementSet())
                weights[i++] = patterns.get(region).count(pattern);
            beagleInstance.setPatternWeights(weights);

            EigenDecomposition ed = substitutionModel.getEigenDecomposition(null);
            beagleInstance.setEigenDecomposition(0,
                    ed.getEigenVectors(),
                    ed.getInverseEigenVectors(),
                    ed.getEigenValues());

            beagleInstance.setCategoryRates(siteModel.getCategoryRates(null));
            beagleInstance.setCategoryWeights(0, siteModel.getCategoryProportions(null));
            beagleInstance.setStateFrequencies(0, substitutionModel.getFrequencies());

            beagleInstances.put(region, beagleInstance);
        }
    }
    
    
    /**
     * Set leaf states in a Beagle instance
     * 
     * @param beagle beagle instance object
     * @param patterns leaf state patterns
     */
    void setStates(Beagle beagle, Multiset<int[]> patterns) {
        
        for (Node node : acg.getExternalNodes()) {
            int[] states = new int[patterns.size()];
            int taxon = alignment.getTaxonIndex(node.getID());
            int i=0;
            for (int [] pattern : patterns.elementSet()) {
                int code = pattern[taxon];
                int[] statesForCode = alignment.getDataType().getStatesForCode(code);
                if (statesForCode.length==1)
                    states[i] = statesForCode[0];
                else
                    states[i] = code; // Causes ambiguous states to be ignored.
                
                i += 1;
            }
            beagle.setTipStates(node.getNr(), states);
        }
    }


    /**
     * Set leaf partials in a Beagle instance
     *
     * @param beagle beagle instance object
     * @param patterns leaf state patterns
     */
    protected void setPartials(Beagle beagle, Multiset<int[]> patterns) {
        for (Node node : acg.getExternalNodes()) {
            Alignment data = dataInput.get();
            int nStates = data.getDataType().getStateCount();
            double[] partials = new double[patterns.elementSet().size() * nStates * siteModel.getCategoryCount()];
            int k = 0;
            int iTaxon = alignment.getTaxonIndex(node.getID());
            for (int[] pattern : patterns.elementSet()) {
                int code = pattern[iTaxon];
                boolean[] stateSet = alignment.getDataType().getStateSet(code);
                for (int iState = 0; iState < nStates; iState++) {
                    partials[k++] = (stateSet[iState] ? 1.0 : 0.0);
                }
            }

            int n = patterns.elementSet().size()*siteModel.getCategoryCount();
            for (int cIdx = 1; cIdx<siteModel.getCategoryCount(); cIdx++) {
                System.arraycopy(partials, 0, partials, n*cIdx, n);
            }

            beagle.setTipPartials(node.getNr(), partials);
        }
    }

    protected void traverse(Beagle beagle, MarginalNode node, Region region) {
        if (!node.isRoot()) {
            edgeLengths[node.getNr()] = node.getLength() * branchRateModel.getRateForBranch(node);
        }

        if (!node.isLeaf()) {

            MarginalNode leftChild = (MarginalNode)node.getLeft();
            MarginalNode rightChild = (MarginalNode)node.getRight();

            traverse(beagle, leftChild, region);
            traverse(beagle, rightChild, region);

            int opIdx = operationListIdx[0];

            operationList[opIdx + 0] = node.getNr();
            operationList[opIdx + 1] = Beagle.NONE;
            operationList[opIdx + 2] = Beagle.NONE;
            operationList[opIdx + 3] = leftChild.getNr();
            operationList[opIdx + 4] = leftChild.getNr();
            operationList[opIdx + 5] = rightChild.getNr();
            operationList[opIdx + 6] = rightChild.getNr();

            operationListIdx[0] += Beagle.OPERATION_TUPLE_SIZE;
        }
    }

    @Override
    public List<String> getArguments() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public List<String> getConditions() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    protected boolean requiresRecalculation() {
        if (acg.clonalFrameIsDirty() || siteModel.isDirtyCalculation())
            regionLogLikelihoods.clear();

        return true;
    }

    @Override
    public void store() {
        storedPatterns.clear();
        storedPatterns.putAll(patterns);

        storedConstantPatterns.clear();
        storedConstantPatterns.putAll(constantPatterns);

        storedRegionLogLikelihoods.clear();
        storedRegionLogLikelihoods.putAll(regionLogLikelihoods);

        super.store();
    }

    @Override
    public void restore() {
        Map<Region, Multiset<int[]>> tmpPatterns = patterns;
        patterns = storedPatterns;
        storedPatterns = tmpPatterns;

        Map<Region, List<Integer>> tmpConstantPatterns = constantPatterns;
        constantPatterns = storedConstantPatterns;
        storedConstantPatterns = tmpConstantPatterns;

        Map<Region, Double> tmpRegionLogLikelihoods = regionLogLikelihoods;
        regionLogLikelihoods = storedRegionLogLikelihoods;
        storedRegionLogLikelihoods = tmpRegionLogLikelihoods;

        super.restore();
    }
}
