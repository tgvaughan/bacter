/*
 * Copyright (C) 2013 Tim Vaughan <tgvaughan@gmail.com>
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
import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.BeerLikelihoodCore;
import beast.evolution.likelihood.BeerLikelihoodCore4;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.likelihood.LikelihoodCore;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import com.google.common.collect.LinkedHashMultiset;
import com.google.common.collect.Multiset;

import java.util.*;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Probability of sequence data given recombination graph.")
public class ACGLikelihood extends GenericTreeLikelihood {

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
    protected BranchRateModel.Base branchRateModel;
    protected SubstitutionModel.Base substitutionModel;
    protected Alignment alignment;
    protected Locus locus;
    protected int nStates;

    protected Map<Region, Multiset<int[]>> patterns;
    protected Map<Region, Multiset<int[]>> storedPatterns;
    protected Map<Region, double[]> patternLogLikelihoods;
    protected Map<Region, double[]> storedPatternLogLikelihoods;
    protected Map<Region, double[]> rootPartials;
    protected Map<Region, double[]> storedRootPartials;
    protected Map<Region, List<Integer>> constantPatterns;
    protected Map<Region, List<Integer>> storedConstantPatterns;
    protected Map<Region, LikelihoodCore> likelihoodCores;
    protected Map<Region, LikelihoodCore> storedLikelihoodCores;
    protected Map<Region, Double> regionLogLikelihoods;
    protected Map<Region, Double> storedRegionLogLikelihoods;

    /**
     * Memory for transition probabilities.
     */
    protected double [] probabilities;

    public ACGLikelihood() {
        // We allow alignments to be specified using Locus objects.
        dataInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {

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
        patternLogLikelihoods = new HashMap<>();
        storedPatternLogLikelihoods = new HashMap<>();
        rootPartials = new HashMap<>();
        storedRootPartials = new HashMap<>();
        constantPatterns = new HashMap<>();
        storedConstantPatterns = new HashMap<>();
        likelihoodCores = new HashMap<>();
        storedLikelihoodCores = new HashMap<>();
        regionLogLikelihoods = new HashMap<>();
        storedRegionLogLikelihoods = new HashMap<>();

        // Allocate transition probability memory:
        // (Only the first nStates*nStates elements are usually used.)
        probabilities = new double[(nStates+1)*(nStates+1)];
    }

    protected double scaleFactor = 1.0;
    protected double scaleFactorMultiplier = 1.01;
    protected double maxScaleFactor = 100.0;

    @Override
    public double calculateLogP() {

        while (true) {
            doLogPCalculation();

            if (logP>Double.NEGATIVE_INFINITY || scaleFactor>=maxScaleFactor)
                return logP;

            scaleFactor *= scaleFactorMultiplier;
            Log.warning.println("Turning on scaling to prevent numeric instability. Scale factor: " + scaleFactor);

            for (LikelihoodCore core : likelihoodCores.values())
                core.setUseScaling(scaleFactor);

            requiresRecalculation();
        }
    }



    protected void doLogPCalculation() {
        updatePatterns();
        updateCores();

        preComputeCFTransitionProbs();

        logP = 0.0;

        regionLogLikelihoods.keySet().retainAll(acg.getRegions(locus));

        for (Region region : acg.getRegions(locus)) {

            if (!regionLogLikelihoods.containsKey(region)) {
                traverseNoRecurse(new MarginalTree(acg, region.activeConversions).getRoot(), region);


                double regionLogP = 0.0;
                int i = 0;
                for (int[] pattern : patterns.get(region).elementSet()) {
                    regionLogP += patternLogLikelihoods.get(region)[i]
                            * patterns.get(region).count(pattern);
                    i += 1;
                }
                regionLogLikelihoods.put(region, regionLogP);

                logP += regionLogP;
            } else {
                logP += regionLogLikelihoods.get(region);
            }
        }

//        System.out.println("Cache hit rate: " + cacheHits/(double)(cacheMisses + cacheHits));
    }

    /**
     * Ensure pattern counts are up to date.
     */
    private void updatePatterns() {
        List<Region> regionList = acg.getRegions(locus);

        // Remove stale pattern sets
        patterns.keySet().retainAll(regionList);
        patternLogLikelihoods.keySet().retainAll(regionList);
        rootPartials.keySet().retainAll(regionList);
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

            // Allocate memory for corresponding log likelihoods and root partials
            patternLogLikelihoods.put(region, new double[patSet.elementSet().size()]);
            rootPartials.put(region, new double[patSet.elementSet().size()*nStates]);

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
     * Initialize likelihood cores.
     */
    private void updateCores() {

        List<Region> regionList = acg.getRegions(locus);
        likelihoodCores.keySet().retainAll(regionList);

        for (Region region : regionList) {

            if (likelihoodCores.containsKey(region))
                continue;

            LikelihoodCore likelihoodCore;
            if (nStates==4)
                likelihoodCore = new BeerLikelihoodCore4();
            else
                likelihoodCore = new BeerLikelihoodCore(nStates);

            likelihoodCores.put(region, likelihoodCore);

            likelihoodCore.initialize(acg.getNodeCount(),
                    patterns.get(region).elementSet().size(),
                    siteModel.getCategoryCount(),
                    true, useAmbiguitiesInput.get());

            if (scaleFactor>1.0)
                likelihoodCore.setUseScaling(scaleFactor);

            if (useAmbiguitiesInput.get())
                setPartials(likelihoodCore, patterns.get(region));
            else
                setStates(likelihoodCore, patterns.get(region));
            
            int intNodeCount = acg.getNodeCount()/2;
            for (int i=0; i<intNodeCount; i++)
                likelihoodCore.createNodePartials(intNodeCount+1+i);
        }
    }
    
    
    /**
     * Set leaf states in a likelihood core.
     * 
     * @param lhc       likelihood core object
     * @param patterns  leaf state patterns
     */
    void setStates(LikelihoodCore lhc, Multiset<int[]> patterns) {
        
        for (Node node : acg.getExternalNodes()) {
            int[] states = new int[patterns.elementSet().size()];
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
            lhc.setNodeStates(node.getNr(), states);
        }
    }


    /**
     * Set leaf partials in likelihood core.
     *
     * @param lhc likelihood core object
     * @param patterns leaf state patterns
     */
    protected void setPartials(LikelihoodCore lhc, Multiset<int[]> patterns) {
        for (Node node : acg.getExternalNodes()) {
            int nStates = alignment.getDataType().getStateCount();
            double[] partials = new double[patterns.elementSet().size() * nStates];
            int k = 0;
            int iTaxon = alignment.getTaxonIndex(node.getID());
            for (int[] pattern : patterns.elementSet()) {
                int code = pattern[iTaxon];
                boolean[] stateSet = alignment.getDataType().getStateSet(code);
                for (int iState = 0; iState < nStates; iState++) {
                    partials[k++] = (stateSet[iState] ? 1.0 : 0.0);
                }
            }
            lhc.setNodePartials(node.getNr(), partials);
        }
    }

    /**
     * Cached transition probabilities for CF edges.
     */
    double [][][] cfTransitionProbs;
    int cacheHits = 0;
    int cacheMisses = 0;

    /**
     * Pre-compute transition probabilities for CF edges.
     */
    void preComputeCFTransitionProbs() {
        if (cfTransitionProbs == null)
            cfTransitionProbs = new double[acg.getNodeCount()-1][siteModel.getCategoryCount()][(nStates+1)*(nStates+1)];

        for (int ni=0; ni<acg.getNodeCount()-1; ni++) {
            Node node = acg.getNode(ni);
            for (int ci=0; ci<siteModel.getCategoryCount(); ci++) {
                double jointBranchRate = siteModel.getRateForCategory(ci, node)
                        * branchRateModel.getRateForBranch(node);
                double parentHeight = node.getParent().getHeight();
                double nodeHeight = node.getHeight();

                substitutionModel.getTransitionProbabilities(
                        node,
                        parentHeight,
                        nodeHeight,
                        jointBranchRate,
                        cfTransitionProbs[ni][ci]);
            }
        }
    }

    Deque<MarginalNode> stack = new ArrayDeque<>();
    MarginalNode[] postOrderNodes;

    /**
     * Assemble list of marginal tree nodes for post-order traversal.
     *
     * @param root root of marginal tree
     */
    void computePostOrder(MarginalNode root) {

        if (postOrderNodes == null)
            postOrderNodes = new MarginalNode[acg.getNodeCount()];

        stack.clear();
        int i = 0;

        stack.push(root);
        while (!stack.isEmpty()) {
            MarginalNode n = stack.pop();
            if (!n.isLeaf()) {
                stack.push((MarginalNode)n.getLeft());
                stack.push((MarginalNode)n.getRight());
            }
            postOrderNodes[postOrderNodes.length-1-(i++)] = n;
        }
    }

    /**
     * Traverse a marginal tree, computing partial likelihoods on the way.
     * This version avoids potentially-expensive recursive function calls.
     *
     * @param root Tree node
     * @param region region
     */
    void traverseNoRecurse(MarginalNode root, Region region) {

        computePostOrder(root);

        LikelihoodCore lhc = likelihoodCores.get(region);

        for (MarginalNode node : postOrderNodes) {

            if (!node.isRoot()) {
                lhc.setNodeMatrixForUpdate(node.getNr());

                boolean cfEdge = node.cfNodeNr>=0
                        && !acg.getNode(node.cfNodeNr).isRoot()
                        && acg.getNode(node.cfNodeNr).getParent().getNr()
                           == ((MarginalNode)node.getParent()).cfNodeNr;

                if (!cfEdge) {
                    cacheMisses += 1;

                    for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                        double jointBranchRate = siteModel.getRateForCategory(i, node)
                                * branchRateModel.getRateForBranch(node);
                        double parentHeight = node.getParent().getHeight();
                        double nodeHeight = node.getHeight();

                        substitutionModel.getTransitionProbabilities(
                                node,
                                parentHeight,
                                nodeHeight,
                                jointBranchRate,
                                probabilities);
                        lhc.setNodeMatrix(node.getNr(), i, probabilities);
                    }
                } else {
                    cacheHits += 1;

                    for (int i=0; i<siteModel.getCategoryCount(); i++) {
                        lhc.setNodeMatrix(node.getNr(), i, cfTransitionProbs[node.cfNodeNr][i]);
                    }
                }
            }

            if (!node.isLeaf()) {

                // LikelihoodCore only supports binary trees.
                List<Node> children = node.getChildren();
                lhc.setNodePartialsForUpdate(node.getNr());
                lhc.setNodeStatesForUpdate(node.getNr());
                lhc.calculatePartials(children.get(0).getNr(),
                        children.get(1).getNr(), node.getNr());

                if (node.isRoot()) {
                    double[] frequencies = substitutionModel.getFrequencies();
                    double[] proportions = siteModel.getCategoryProportions(node);
                    lhc.integratePartials(node.getNr(), proportions,
                            rootPartials.get(region));

                    for (int idx : constantPatterns.get(region)) {
                        rootPartials.get(region)[idx]
                                += siteModel.getProportionInvariant();
                    }

                    lhc.calculateLogLikelihoods(rootPartials.get(region),
                            frequencies, patternLogLikelihoods.get(region));
                }
            }
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
        if (acg.clonalFrameIsDirty()
                || siteModel.isDirtyCalculation()
                || branchRateModel.isDirtyCalculation())
            regionLogLikelihoods.clear();

        return true;
    }

    @Override
    public void store() {
        storedPatterns.clear();
        storedPatterns.putAll(patterns);

        storedPatternLogLikelihoods.clear();
        storedPatternLogLikelihoods.putAll(patternLogLikelihoods);

        storedConstantPatterns.clear();
        storedConstantPatterns.putAll(constantPatterns);

        storedRootPartials.clear();
        storedRootPartials.putAll(rootPartials);

        storedLikelihoodCores.clear();
        storedLikelihoodCores.putAll(likelihoodCores);

        storedRegionLogLikelihoods.clear();
        storedRegionLogLikelihoods.putAll(regionLogLikelihoods);

        super.store();
    }

    @Override
    public void restore() {
        Map<Region, Multiset<int[]>> tmpPatterns = patterns;
        patterns = storedPatterns;
        storedPatterns = tmpPatterns;

        Map<Region, double[]> tmpPatternLogLikelihoods = patternLogLikelihoods;
        patternLogLikelihoods = storedPatternLogLikelihoods;
        storedPatternLogLikelihoods = tmpPatternLogLikelihoods;

        Map<Region, double[]> tmpRootPartials = rootPartials;
        rootPartials = storedRootPartials;
        storedRootPartials = tmpRootPartials;

        Map<Region, LikelihoodCore> tmpLikelihoodCores = likelihoodCores;
        likelihoodCores = storedLikelihoodCores;
        storedLikelihoodCores = tmpLikelihoodCores;

        Map<Region, List<Integer>> tmpConstantPatterns = constantPatterns;
        constantPatterns = storedConstantPatterns;
        storedConstantPatterns = tmpConstantPatterns;

        Map<Region, Double> tmpRegionLogLikelihoods = regionLogLikelihoods;
        regionLogLikelihoods = storedRegionLogLikelihoods;
        storedRegionLogLikelihoods = tmpRegionLogLikelihoods;

        super.restore();
    }
}
