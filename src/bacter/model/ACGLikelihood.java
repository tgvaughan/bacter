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

import bacter.Conversion;
import bacter.Locus;
import bacter.MarginalTree;
import bacter.Region;
import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.BeerLikelihoodCore;
import beast.evolution.likelihood.BeerLikelihoodCore4;
import beast.evolution.likelihood.LikelihoodCore;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import com.google.common.collect.LinkedHashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Multiset;

import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Probability of sequence data given recombination graph.")
public class ACGLikelihood extends ACGDistribution {

    public Input<Locus> locusInput = new Input<>(
            "locus",
            "Locus associated with alignment to evaluate probability of.",
            Input.Validate.REQUIRED);

    public Input<Alignment> alignmentInput = new Input<>(
            "alignment",
            "Used to explicitly specify alignment when none is " +
                    "directly associated with given locus");

    public Input<SiteModelInterface> siteModelInput = new Input<>(
            "siteModel",
            "Site model for evolution of alignment.",
            Input.Validate.REQUIRED);

    protected SiteModel.Base siteModel;
    protected SubstitutionModel.Base substitutionModel;
    protected Alignment alignment;
    protected Locus locus;
    protected int nStates;

    protected Map<Set<Conversion>, LikelihoodCore> likelihoodCores;
    
    protected Map<Set<Conversion>, Multiset<int[]>> patterns;
    protected Map<Set<Conversion>, double[]> patternLogLikelihoods;
    protected Map<Set<Conversion>, double[]> rootPartials;
    protected Map<Set<Conversion>, List<Integer>> constantPatterns;
    
    ExecutorService likelihoodThreadPool;
    int threadPoolSize;

    /**
     * Memory for transition probabilities.
     */
    double [] probabilities;
    
    @Override
    public void initAndValidate() throws Exception {

        super.initAndValidate();
        
        acg = acgInput.get();

        locus = locusInput.get();
        if (locus.hasAlignment()) {
            alignment = locus.getAlignment();
        } else {
            if (alignmentInput.get() != null)
                alignment = alignmentInput.get();
            else
                throw new IllegalArgumentException("No alignment associated with " +
                        "locus " + locus.getID() + " provided to ACGLikelihood " +
                        "and none given explicitly.");
        }

        nStates = alignment.getMaxStateCount();

        siteModel = (SiteModel.Base) siteModelInput.get();
        substitutionModel = (SubstitutionModel.Base) siteModel.getSubstitutionModel();

        // Initialize patterns
        patterns = new LinkedHashMap<>();
        patternLogLikelihoods = new HashMap<>();
        rootPartials = new HashMap<>();
        constantPatterns = new HashMap<>();

        // Initialise cores        
        likelihoodCores = new HashMap<>();

        // Allocate transition probability memory:
        // (Only the first nStates*nStates elements are usually used.)
        probabilities = new double[(nStates+1)*(nStates+1)];

        // Initialize thread pool using number of available processors.
        threadPoolSize = Runtime.getRuntime().availableProcessors();
        likelihoodThreadPool = Executors.newFixedThreadPool(threadPoolSize);
    }
    
    
    @Override
    public double calculateACGLogP() {
        
        logP = 0.0;

        updatePatterns();
        updateCores();

        int nRegions = patterns.size();
        int usedPoolSize = Integer.min(threadPoolSize, nRegions);
        int chunkSize =  nRegions % usedPoolSize == 0
                ? nRegions/usedPoolSize
                : nRegions/usedPoolSize + 1;


        List<Set<Conversion>> convSets = new ArrayList<>(patterns.keySet());
        Future<?>[] futures = new Future[usedPoolSize];

        // Submit likelihood calculations to thread pool members
        for (int threadIdx=0; threadIdx<usedPoolSize; threadIdx++) {

            int start = threadIdx*chunkSize;
            int end = Integer.min((threadIdx+1)*chunkSize, convSets.size());

            if (end<=start)
                continue;

            List<Set<Conversion>> theseConvSets = convSets.subList(start, end);

            futures[threadIdx] = likelihoodThreadPool.submit(() -> {
                for (Set<Conversion> convSet : theseConvSets)
                    traverse(new MarginalTree(acg, convSet).getRoot(), convSet);
            });
        }

        // Wait for all tasks to complete.
        try {
            for (int i=0; i<usedPoolSize; i++)
                if (futures[i] != null)
                    futures[i].get();
        } catch (InterruptedException | ExecutionException e) {
            throw new IllegalStateException("ACG likelihood calculation task interrupted.");
        }

        // Collate log likelihoods
        for (Set<Conversion> convSet : patterns.keySet()) {
            int i=0;
            for (int[] pattern : patterns.get(convSet).elementSet()) {
                logP += patternLogLikelihoods.get(convSet)[i]
                        *patterns.get(convSet).count(pattern);
                i += 1;
            }
        }

        return logP;
    }
    
    
    /**
     * Ensure pattern counts are up to date.
     */
    private void updatePatterns() {

        patterns.clear();
        patternLogLikelihoods.clear();
        rootPartials.clear();
        
        for (Region region : acg.getRegions(locus)) {

            Multiset<int[]> patSet;
            if (patterns.containsKey(region.activeConversions))
                patSet = patterns.get(region.activeConversions);
            else {
                patSet = LinkedHashMultiset.create();
                patterns.put(region.activeConversions, patSet);
            }

            for (int j=region.leftBoundary; j<region.rightBoundary; j++) {
                int [] pat = alignment.getPattern(alignment.getPatternIndex(j));
                patSet.add(pat);
            }
        }

        // Allocate memory for likelihoods and partials, and construct list
        // of constant patterns.
        constantPatterns.clear();
        for (Set<Conversion> convSet : patterns.keySet()) {
            Multiset<int[]> patSet = patterns.get(convSet);
            patternLogLikelihoods.put(convSet, new double[patSet.elementSet().size()]);
            rootPartials.put(convSet, new double[patSet.elementSet().size()*nStates]);

            List<Integer> constantPatternList = Lists.newArrayList();
            
            int patternIdx = 0;
            for (int[] pattern : patterns.get(convSet).elementSet()) {
                boolean isConstant = true;
                for (int i=1; i<pattern.length; i++)
                    if (pattern[i] != pattern[0]) {
                        isConstant = false;
                        break;
                    }
                
                if (isConstant && !alignment.getDataType().isAmbiguousState(pattern[0]))
                    constantPatternList.add(patternIdx*nStates + pattern[0]);
                
                patternIdx += 1;
            }
            
            constantPatterns.put(convSet, constantPatternList);
       }
    }
    
    
    /**
     * Initialize likelihood cores.
     */
    private void updateCores() {
        
        likelihoodCores.clear();
        
        for (Set<Conversion> convSet : patterns.keySet()) {
            
            LikelihoodCore likelihoodCore;
            if (!likelihoodCores.keySet().contains(convSet)) {
                if (nStates==4)
                    likelihoodCore = new BeerLikelihoodCore4();
                else
                    likelihoodCore = new BeerLikelihoodCore(nStates);
                
                likelihoodCores.put(convSet, likelihoodCore);
            } else
                likelihoodCore = likelihoodCores.get(convSet);
            
            likelihoodCore.initialize(acg.getNodeCount(),
                patterns.get(convSet).elementSet().size(),
                siteModel.getCategoryCount(),
                true, false);
            setStates(likelihoodCore, patterns.get(convSet));
            
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
            lhc.setNodeStates(node.getNr(), states);
        }
    }
    
    
    /**
     * Traverse a marginal tree, computing partial likelihoods on the way.
     * 
     * @param node Tree node
     * @param convSet Set of active conversions.
     */
    void traverse(Node node, Set<Conversion> convSet) {

        LikelihoodCore lhc = likelihoodCores.get(convSet);
        
        if (!node.isRoot()) {
            lhc.setNodeMatrixForUpdate(node.getNr());
            for (int i=0; i<siteModel.getCategoryCount(); i++) {
                double jointBranchRate = siteModel.getRateForCategory(i, node);
                double parentHeight = node.getParent().getHeight();
                double nodeHeight = node.getHeight();

                synchronized (this) {
                    substitutionModel.getTransitionProbabilities(
                            node,
                            parentHeight,
                            nodeHeight,
                            jointBranchRate,
                            probabilities);
                    lhc.setNodeMatrix(node.getNr(), i, probabilities);
                }
            }
        }
        
        if (!node.isLeaf()) {
            
            // LikelihoodCore only supports binary trees.
            List<Node> children = node.getChildren();
            traverse(children.get(0), convSet);
            traverse(children.get(1), convSet);

            lhc.setNodePartialsForUpdate(node.getNr());
            lhc.setNodeStatesForUpdate(node.getNr());
            lhc.calculatePartials(children.get(0).getNr(),
                    children.get(1).getNr(), node.getNr());

            if (node.isRoot()) {
                double [] frequencies = substitutionModel.getFrequencies();
                double [] proportions = siteModel.getCategoryProportions(node);
                lhc.integratePartials(node.getNr(), proportions,
                        rootPartials.get(convSet));

                for (int idx : constantPatterns.get(convSet)) {
                    rootPartials.get(convSet)[idx]
                            += siteModel.getProportionInvariant();
                }

                lhc.calculateLogLikelihoods(rootPartials.get(convSet),
                        frequencies, patternLogLikelihoods.get(convSet));
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
}
