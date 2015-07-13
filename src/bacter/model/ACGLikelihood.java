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

import bacter.ConversionGraph;
import bacter.Locus;
import bacter.MarginalTree;
import bacter.Region;
import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.evolution.alignment.Alignment;
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
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Probability of sequence data given recombination graph.")
public class ACGLikelihood extends GenericTreeLikelihood {

    public Input<Locus> locusInput = new Input<>(
            "locus",
            "Locus associated with alignment to evaluate probability of.",
            Input.Validate.REQUIRED);

    public Input<Integer> nThreadsInput = new Input<>(
            "nThreads",
            "Number of additional threads to use in calculation.",
            0);

    protected ConversionGraph acg;

    protected SiteModel.Base siteModel;
    protected SubstitutionModel.Base substitutionModel;
    protected Alignment alignment;
    protected Locus locus;
    protected int nStates;

    protected Map<Region, LikelihoodCore> likelihoodCores;
    
    protected Map<Region, Multiset<int[]>> patterns;
    protected Map<Region, double[]> patternLogLikelihoods;
    protected Map<Region, double[]> rootPartials;
    protected Map<Region, List<Integer>> constantPatterns;

    ExecutorService likelihoodThreadPool;

    /**
     * Memory for transition probabilities.
     */
    double [] probabilities;

    public ACGLikelihood() {
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

        if (branchRateModelInput.get() != null &&
                !(branchRateModelInput.get() instanceof StrictClockModel))
            throw new IllegalArgumentException("ACGLikelihood currently only" +
                    "supports strict clock models.");

        // Initialize patterns
//        patterns = new TreeMap<>((r1, r2) -> {
//            if (r1.leftBoundary<r2.leftBoundary)
//                return -1;
//            if (r1.leftBoundary>r1.rightBoundary)
//                return 1;
//            return 0;
//        });
//        patternLogLikelihoods = new TreeMap<>((r1, r2) -> {
//        if (r1.leftBoundary<r2.leftBoundary)
//            return -1;
//            if (r1.leftBoundary>r1.rightBoundary)
//                return 1;
//            return 0;
//        });
//        rootPartials = new TreeMap<>((r1, r2) -> {
//            if (r1.leftBoundary<r2.leftBoundary)
//                return -1;
//            if (r1.leftBoundary>r1.rightBoundary)
//                return 1;
//            return 0;
//        });
//        constantPatterns = new TreeMap<>((r1, r2) -> {
//            if (r1.leftBoundary<r2.leftBoundary)
//                return -1;
//            if (r1.leftBoundary>r1.rightBoundary)
//                return 1;
//            return 0;
//        });
//
//        // Initialise cores
//        likelihoodCores = new TreeMap<>((r1, r2) -> {
//            if (r1.leftBoundary<r2.leftBoundary)
//                return -1;
//            if (r1.leftBoundary>r1.rightBoundary)
//                return 1;
//            return 0;
//        });

        patterns = new HashMap<>();
        patternLogLikelihoods = new HashMap<>();
        rootPartials = new HashMap<>();
        constantPatterns = new HashMap<>();
        likelihoodCores = new HashMap<>();

        // Allocate transition probability memory:
        // (Only the first nStates*nStates elements are usually used.)
        probabilities = new double[(nStates+1)*(nStates+1)];

        // Initialize thread pool using number of available processors.
        if (nThreadsInput.get()>0)
            likelihoodThreadPool = Executors.newFixedThreadPool(nThreadsInput.get());
    }

    @Override
    public double calculateLogP() {
//        updatePatterns();
//        updateCores();
//
//        if (nThreadsInput.get() > 0)
//            calculateLogPThreads();
//        else
//            calculateLogPNoThreads();
//
//        double logP1 = logP;

        //clearPatterns();

        updatePatterns();
        updateCores();

        if (nThreadsInput.get() > 0)
            calculateLogPThreads();
        else
            calculateLogPNoThreads();

        double logP2 = logP;

        return logP;
    }


    /**
     * Compute likelihood using single thread.
     *
     * @return logP
     */
    public double calculateLogPNoThreads() {
        logP = 0.0;

        for (Region region : acg.getRegions(locus)) {
            traverse(new MarginalTree(acg, region.activeConversions).getRoot(), region);

            int i=0;
            for (int[] pattern : patterns.get(region).elementSet()) {
                logP += patternLogLikelihoods.get(region)[i]
                        *patterns.get(region).count(pattern);
                i += 1;
            }
        }

        return logP;
    }

    /**
     * Compute likelihood using multiple threads.  This allows different
     * marginal trees to be considered in parallel.
     *
     * @return logP
     */
    public double calculateLogPThreads() {
        logP = 0.0;

        int nRegions = patterns.size();
        int usedPoolSize = Integer.min(nThreadsInput.get(), nRegions);
        int chunkSize =  nRegions % usedPoolSize == 0
                ? nRegions/usedPoolSize
                : nRegions/usedPoolSize + 1;


        List<Region> regions = new ArrayList<>(patterns.keySet());
        Future<?>[] futures = new Future[usedPoolSize];

        // Submit likelihood calculations to thread pool members
        for (int threadIdx=0; threadIdx<usedPoolSize; threadIdx++) {

            int start = threadIdx*chunkSize;
            int end = Integer.min((threadIdx + 1) * chunkSize, regions.size());

            if (end<=start)
                continue;

            List<Region> theseRegions = regions.subList(start, end);

            futures[threadIdx] = likelihoodThreadPool.submit(() -> {
                for (Region region : theseRegions)
                    traverse(new MarginalTree(acg, region.activeConversions).getRoot(),
                            region);
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
        for (Region region : patterns.keySet()) {
            int i=0;
            for (int[] pattern : patterns.get(region).elementSet()) {
                logP += patternLogLikelihoods.get(region)[i]
                        *patterns.get(region).count(pattern);
                i += 1;
            }
        }

        return logP;
    }

    private void clearPatterns() {
        patterns.clear();
        constantPatterns.clear();
        patternLogLikelihoods.clear();
        rootPartials.clear();
        likelihoodCores.clear();
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
            for (int[] pattern : patterns.get(region).elementSet()) {
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
                    true, false);
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
     * @param region region
     */
    void traverse(Node node, Region region) {

        LikelihoodCore lhc = likelihoodCores.get(region);
        
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
            traverse(children.get(0), region);
            traverse(children.get(1), region);

            lhc.setNodePartialsForUpdate(node.getNr());
            lhc.setNodeStatesForUpdate(node.getNr());
            lhc.calculatePartials(children.get(0).getNr(),
                    children.get(1).getNr(), node.getNr());

            if (node.isRoot()) {
                double [] frequencies = substitutionModel.getFrequencies();
                double [] proportions = siteModel.getCategoryProportions(node);
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
