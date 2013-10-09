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
package argbeast;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.likelihood.BeerLikelihoodCore;
import beast.evolution.likelihood.BeerLikelihoodCore4;
import beast.evolution.likelihood.LikelihoodCore;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.util.ClusterTree;
import com.google.common.collect.LinkedHashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multiset;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Probability of sequence data given recombination graph.")
public class RecombinationGraphLikelihood extends Distribution {
    
    public Input<RecombinationGraph> argInput = new Input<RecombinationGraph>(
            "arg", "Recombination graph.", Validate.REQUIRED);
    
    public Input<Alignment> alignmentInput = new Input<Alignment>(
            "data", "Sequence data to evaluate probability of.",
            Validate.REQUIRED);
    
    public Input<SiteModelInterface> siteModelInput = new Input<SiteModelInterface>(
            "siteModel", "Site model for evolution of alignment.",
            Validate.REQUIRED);

    RecombinationGraph arg;
    SiteModel.Base siteModel;
    SubstitutionModel.Base substitutionModel;
    Alignment alignment;
    
    Map<Recombination,LikelihoodCore> likelihoodCores;
    
    Map<Recombination, Multiset<int[]>> patterns;
    Map<Recombination, double[]> patternLogLikelihoods;
    Map<Recombination, double[]> rootPartials;
    Map<Recombination, List<Integer>> constantPatterns;
    
    int nStates;
    
    /*
     * Memory for transition probabilities.
     */
    double [] probabilities;
    
    @Override
    public void initAndValidate() throws Exception {
        
        arg = argInput.get();
        alignment = alignmentInput.get();
        
        nStates = alignment.getMaxStateCount();

        // Initialize patterns
        patterns = Maps.newHashMap();
        updatePatterns();
        
        // Initialise cores        
        likelihoodCores = Maps.newHashMap();
        updateCores();
        
        // Allocate transition probability memory:
        // (Only the first nStates*nStates elements are usually used.)
        probabilities = new double[(nStates+1)*(nStates+1)];
    }
    
    
    /**
     * Ensure pattern counts are up to date.
     */
    private void updatePatterns() {

        patterns.clear();
        patternLogLikelihoods.clear();
        rootPartials.clear();
        
        Multiset<int[]> cfPatSet = LinkedHashMultiset.create();

        int j=0;
        for (Recombination recomb : arg.getRecombinations()) {
            if (recomb == null)
                continue; // Skip clonal frame
            
            while (j < recomb.startLocus) {
                int [] pat = alignment.getPattern(alignment.getPatternIndex(j));
                cfPatSet.add(pat);
                j += 1;
            }
            
            Multiset<int[]> recombPatSet = LinkedHashMultiset.create();
            
            while (j <= recomb.endLocus) {                
                int [] pat = alignment.getPattern(alignment.getPatternIndex(j));
                recombPatSet.add(pat);
                j += 1;
            }
            
            patterns.put(recomb, recombPatSet);
            patternLogLikelihoods.put(recomb,
                    new double[recombPatSet.elementSet().size()]);
            rootPartials.put(recomb,
                    new double[recombPatSet.elementSet().size()*nStates]);
        }
        
        while (j<alignmentInput.get().getSiteCount()) {
            
            int [] pat = alignment.getPattern(alignment.getPatternIndex(j));
            cfPatSet.add(pat);
            j += 1;
        }

        patterns.put(null, cfPatSet);
        patternLogLikelihoods.put(null,
                new double[cfPatSet.elementSet().size()]);
        rootPartials.put(null,
                    new double[cfPatSet.elementSet().size()*nStates]);
        
        // Record lists of constant patterns:
        constantPatterns.clear();
        for (Recombination recomb: arg.getRecombinations()) {
            List<Integer> constantPatternList = Lists.newArrayList();
            
            int patternIdx = 0;
            for (int[] pattern : patterns.get(recomb).elementSet()) {
                boolean isConstant = true;
                for (int i=1; i<pattern.length; i++)
                    if (pattern[i] != pattern[0]) {
                        isConstant = false;
                        break;
                    }
                
                if (isConstant)
                    constantPatternList.add(patternIdx + pattern[0]);
                
                patternIdx += 1;
            }
            
            constantPatterns.put(recomb, constantPatternList);
        }
    }
    
    
    /**
     * Initialise likelihood cores.
     */
    private void updateCores() {
        
        likelihoodCores.keySet().retainAll(arg.recombs);
        
        for (Recombination recomb : arg.getRecombinations()) {
            
            LikelihoodCore likelihoodCore;
            if (!likelihoodCores.keySet().contains(recomb)) {
                if (nStates==4)
                    likelihoodCore = new BeerLikelihoodCore4();
                else
                    likelihoodCore = new BeerLikelihoodCore(nStates);
                
                likelihoodCores.put(recomb, likelihoodCore);
            } else
                likelihoodCore = likelihoodCores.get(recomb);
            
            likelihoodCore.initialize(
                arg.getNodeCount(),
                patterns.get(recomb).elementSet().size(),
                siteModel.getCategoryCount(),
                true, false);
            setStates(likelihoodCore, patterns.get(recomb));
        }
    }
    
    
    /**
     * Set leaf states in a likelihood core.
     * 
     * @param lhc
     * @param patterns 
     */
    void setStates(LikelihoodCore lhc, Multiset<int[]> patterns) {
        
        for (Node node : arg.getExternalNodes()) {
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
    
    @Override
    public double calculateLogP() {
        
        logP = 0.0;
        
        updatePatterns();
        updateCores();
        
        for (Recombination recomb : arg.getRecombinations()) {
            traverse(arg.getMarginalRoot(recomb), recomb);
            
            int i=0;
            for (int[] pattern : patterns.get(recomb)) {
                logP += patternLogLikelihoods.get(recomb)[i]
                        *patterns.get(recomb).count(pattern);
                i += 1;
            }
        }
        
        return logP;
    }
    
    
    /**
     * Traverse a marginal tree, computing partial likelihoods on the way.
     * 
     * @param node Tree node
     * @param recomb Recombination object.  Null selects the clonal frame.
     */
    void traverse(Node node, Recombination recomb) {

        LikelihoodCore lhc = likelihoodCores.get(recomb);
        
        double branchTime = arg.getMarginalBranchLength(node, recomb);
        
        if (!arg.isNodeMarginalRoot(node, recomb)) {
            lhc.setNodeMatrixForUpdate(node.getNr());
            for (int i=0; i<siteModel.getCategoryCount(); i++) {
                double jointBranchRate = siteModel.getRateForCategory(i, node);
                double parentHeight = arg.getMarginalNodeHeight(
                        arg.getMarginalParent(node, recomb), recomb);
                double nodeHeight = arg.getMarginalNodeHeight(node, recomb);
                substitutionModel.getTransitionProbabilities(
                        node,
                        parentHeight,
                        nodeHeight,
                        jointBranchRate,
                        probabilities);
            }
        }
        
        if (!arg.isNodeMarginalLeaf(node, recomb)) {
            
            // LikelihoodCore only supports binary trees.
            List<Node> children = arg.getMarginalChildren(node, recomb);
            traverse(children.get(0), recomb);
            traverse(children.get(1), recomb);
            
            lhc.setNodePartialsForUpdate(node.getNr());
            lhc.setNodeStatesForUpdate(node.getNr());
            lhc.calculatePartials(children.get(0).getNr(),
                    children.get(1).getNr(), node.getNr());
            
            if (arg.isNodeMarginalRoot(node, recomb)) {
                double [] frequencies = substitutionModel.getFrequencies();
                double [] proportions = siteModel.getCategoryProportions(node);
                lhc.integratePartials(node.getNr(), proportions,
                        rootPartials.get(recomb));
                
                for (int idx : constantPatterns.get(recomb)) {
                    rootPartials.get(recomb)[idx]
                            += siteModel.getProportionInvariant();
                }
                
                lhc.calculateLogLikelihoods(rootPartials.get(recomb),
                        frequencies, patternLogLikelihoods.get(recomb));
            }
        }
    }
    
    
    @Override
    public List<String> getArguments() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public List<String> getConditions() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    
    /**
     * Main method for testing.
     * 
     * @param args 
     */
    public static void main(String[] args) throws Exception {
       
        // Sequence alignment
        
        List<Sequence> sequences = new ArrayList<Sequence>();
        sequences.add(new Sequence("Tarsius_syrichta","AAGTTTCATTGGAGCCACCACTCTTATAATTGCCCATGGCCTCACCTCCTCCCTATTATTTTGCCTAGCAAATACAAACTACGAACGAGTCCACAGTCGAACAATAGCACTAGCCCGTGGCCTTCAAACCCTATTACCTCTTGCAGCAACATGATGACTCCTCGCCAGCTTAACCAACCTGGCCCTTCCCCCAACAATTAATTTAATCGGTGAACTGTCCGTAATAATAGCAGCATTTTCATGGTCACACCTAACTATTATCTTAGTAGGCCTTAACACCCTTATCACCGCCCTATATTCCCTATATATACTAATCATAACTCAACGAGGAAAATACACATATCATATCAACAATATCATGCCCCCTTTCACCCGAGAAAATACATTAATAATCATACACCTATTTCCCTTAATCCTACTATCTACCAACCCCAAAGTAATTATAGGAACCATGTACTGTAAATATAGTTTAAACAAAACATTAGATTGTGAGTCTAATAATAGAAGCCCAAAGATTTCTTATTTACCAAGAAAGTA-TGCAAGAACTGCTAACTCATGCCTCCATATATAACAATGTGGCTTTCTT-ACTTTTAAAGGATAGAAGTAATCCATCGGTCTTAGGAACCGAAAA-ATTGGTGCAACTCCAAATAAAAGTAATAAATTTATTTTCATCCTCCATTTTACTATCACTTACACTCTTAATTACCCCATTTATTATTACAACAACTAAAAAATATGAAACACATGCATACCCTTACTACGTAAAAAACTCTATCGCCTGCGCATTTATAACAAGCCTAGTCCCAATGCTCATATTTCTATACACAAATCAAGAAATAATCATTTCCAACTGACATTGAATAACGATTCATACTATCAAATTATGCCTAAGCTT"));
        sequences.add(new Sequence("Lemur_catta","AAGCTTCATAGGAGCAACCATTCTAATAATCGCACATGGCCTTACATCATCCATATTATTCTGTCTAGCCAACTCTAACTACGAACGAATCCATAGCCGTACAATACTACTAGCACGAGGGATCCAAACCATTCTCCCTCTTATAGCCACCTGATGACTACTCGCCAGCCTAACTAACCTAGCCCTACCCACCTCTATCAATTTAATTGGCGAACTATTCGTCACTATAGCATCCTTCTCATGATCAAACATTACAATTATCTTAATAGGCTTAAATATGCTCATCACCGCTCTCTATTCCCTCTATATATTAACTACTACACAACGAGGAAAACTCACATATCATTCGCACAACCTAAACCCATCCTTTACACGAGAAAACACCCTTATATCCATACACATACTCCCCCTTCTCCTATTTACCTTAAACCCCAAAATTATTCTAGGACCCACGTACTGTAAATATAGTTTAAA-AAAACACTAGATTGTGAATCCAGAAATAGAAGCTCAAAC-CTTCTTATTTACCGAGAAAGTAATGTATGAACTGCTAACTCTGCACTCCGTATATAAAAATACGGCTATCTCAACTTTTAAAGGATAGAAGTAATCCATTGGCCTTAGGAGCCAAAAA-ATTGGTGCAACTCCAAATAAAAGTAATAAATCTATTATCCTCTTTCACCCTTGTCACACTGATTATCCTAACTTTACCTATCATTATAAACGTTACAAACATATACAAAAACTACCCCTATGCACCATACGTAAAATCTTCTATTGCATGTGCCTTCATCACTAGCCTCATCCCAACTATATTATTTATCTCCTCAGGACAAGAAACAATCATTTCCAACTGACATTGAATAACAATCCAAACCCTAAAACTATCTATTAGCTT"));
        sequences.add(new Sequence("Homo_sapiens","AAGCTTCACCGGCGCAGTCATTCTCATAATCGCCCACGGGCTTACATCCTCATTACTATTCTGCCTAGCAAACTCAAACTACGAACGCACTCACAGTCGCATCATAATCCTCTCTCAAGGACTTCAAACTCTACTCCCACTAATAGCTTTTTGATGACTTCTAGCAAGCCTCGCTAACCTCGCCTTACCCCCCACTATTAACCTACTGGGAGAACTCTCTGTGCTAGTAACCACGTTCTCCTGATCAAATATCACTCTCCTACTTACAGGACTCAACATACTAGTCACAGCCCTATACTCCCTCTACATATTTACCACAACACAATGGGGCTCACTCACCCACCACATTAACAACATAAAACCCTCATTCACACGAGAAAACACCCTCATGTTCATACACCTATCCCCCATTCTCCTCCTATCCCTCAACCCCGACATCATTACCGGGTTTTCCTCTTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTTA-CGACCCCTTATTTACCGAGAAAGCT-CACAAGAACTGCTAACTCATGCCCCCATGTCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGCACACTACTATAACCACCCTAACCCTGACTTCCCTAATTCCCCCCATCCTTACCACCCTCGTTAACCCTAACAAAAAAAACTCATACCCCCATTATGTAAAATCCATTGTCGCATCCACCTTTATTATCAGTCTCTTCCCCACAACAATATTCATGTGCCTAGACCAAGAAGTTATTATCTCGAACTGACACTGAGCCACAACCCAAACAACCCAGCTCTCCCTAAGCTT"));
        sequences.add(new Sequence("Pan","AAGCTTCACCGGCGCAATTATCCTCATAATCGCCCACGGACTTACATCCTCATTATTATTCTGCCTAGCAAACTCAAATTATGAACGCACCCACAGTCGCATCATAATTCTCTCCCAAGGACTTCAAACTCTACTCCCACTAATAGCCTTTTGATGACTCCTAGCAAGCCTCGCTAACCTCGCCCTACCCCCTACCATTAATCTCCTAGGGGAACTCTCCGTGCTAGTAACCTCATTCTCCTGATCAAATACCACTCTCCTACTCACAGGATTCAACATACTAATCACAGCCCTGTACTCCCTCTACATGTTTACCACAACACAATGAGGCTCACTCACCCACCACATTAATAACATAAAGCCCTCATTCACACGAGAAAATACTCTCATATTTTTACACCTATCCCCCATCCTCCTTCTATCCCTCAATCCTGATATCATCACTGGATTCACCTCCTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTCA-CGACCCCTTATTTACCGAGAAAGCT-TATAAGAACTGCTAATTCATATCCCCATGCCTGACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCCATCCGTTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGTATACTACCATAACCACCTTAACCCTAACTCCCTTAATTCTCCCCATCCTCACCACCCTCATTAACCCTAACAAAAAAAACTCATATCCCCATTATGTGAAATCCATTATCGCGTCCACCTTTATCATTAGCCTTTTCCCCACAACAATATTCATATGCCTAGACCAAGAAGCTATTATCTCAAACTGGCACTGAGCAACAACCCAAACAACCCAGCTCTCCCTAAGCTT"));
        sequences.add(new Sequence("Gorilla","AAGCTTCACCGGCGCAGTTGTTCTTATAATTGCCCACGGACTTACATCATCATTATTATTCTGCCTAGCAAACTCAAACTACGAACGAACCCACAGCCGCATCATAATTCTCTCTCAAGGACTCCAAACCCTACTCCCACTAATAGCCCTTTGATGACTTCTGGCAAGCCTCGCCAACCTCGCCTTACCCCCCACCATTAACCTACTAGGAGAGCTCTCCGTACTAGTAACCACATTCTCCTGATCAAACACCACCCTTTTACTTACAGGATCTAACATACTAATTACAGCCCTGTACTCCCTTTATATATTTACCACAACACAATGAGGCCCACTCACACACCACATCACCAACATAAAACCCTCATTTACACGAGAAAACATCCTCATATTCATGCACCTATCCCCCATCCTCCTCCTATCCCTCAACCCCGATATTATCACCGGGTTCACCTCCTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGATAACAGAGGCTCA-CAACCCCTTATTTACCGAGAAAGCT-CGTAAGAGCTGCTAACTCATACCCCCGTGCTTGACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACTATGTACGCTACCATAACCACCTTAGCCCTAACTTCCTTAATTCCCCCTATCCTTACCACCTTCATCAATCCTAACAAAAAAAGCTCATACCCCCATTACGTAAAATCTATCGTCGCATCCACCTTTATCATCAGCCTCTTCCCCACAACAATATTTCTATGCCTAGACCAAGAAGCTATTATCTCAAGCTGACACTGAGCAACAACCCAAACAATTCAACTCTCCCTAAGCTT"));
        sequences.add(new Sequence("Pongo","AAGCTTCACCGGCGCAACCACCCTCATGATTGCCCATGGACTCACATCCTCCCTACTGTTCTGCCTAGCAAACTCAAACTACGAACGAACCCACAGCCGCATCATAATCCTCTCTCAAGGCCTTCAAACTCTACTCCCCCTAATAGCCCTCTGATGACTTCTAGCAAGCCTCACTAACCTTGCCCTACCACCCACCATCAACCTTCTAGGAGAACTCTCCGTACTAATAGCCATATTCTCTTGATCTAACATCACCATCCTACTAACAGGACTCAACATACTAATCACAACCCTATACTCTCTCTATATATTCACCACAACACAACGAGGTACACCCACACACCACATCAACAACATAAAACCTTCTTTCACACGCGAAAATACCCTCATGCTCATACACCTATCCCCCATCCTCCTCTTATCCCTCAACCCCAGCATCATCGCTGGGTTCGCCTACTGTAAATATAGTTTAACCAAAACATTAGATTGTGAATCTAATAATAGGGCCCCA-CAACCCCTTATTTACCGAGAAAGCT-CACAAGAACTGCTAACTCTCACT-CCATGTGTGACAACATGGCTTTCTCAGCTTTTAAAGGATAACAGCTATCCCTTGGTCTTAGGATCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAACAGCCATGTTTACCACCATAACTGCCCTCACCTTAACTTCCCTAATCCCCCCCATTACCGCTACCCTCATTAACCCCAACAAAAAAAACCCATACCCCCACTATGTAAAAACGGCCATCGCATCCGCCTTTACTATCAGCCTTATCCCAACAACAATATTTATCTGCCTAGGACAAGAAACCATCGTCACAAACTGATGCTGAACAACCACCCAGACACTACAACTCTCACTAAGCTT"));
        sequences.add(new Sequence("Hylobates","AAGCTTTACAGGTGCAACCGTCCTCATAATCGCCCACGGACTAACCTCTTCCCTGCTATTCTGCCTTGCAAACTCAAACTACGAACGAACTCACAGCCGCATCATAATCCTATCTCGAGGGCTCCAAGCCTTACTCCCACTGATAGCCTTCTGATGACTCGCAGCAAGCCTCGCTAACCTCGCCCTACCCCCCACTATTAACCTCCTAGGTGAACTCTTCGTACTAATGGCCTCCTTCTCCTGGGCAAACACTACTATTACACTCACCGGGCTCAACGTACTAATCACGGCCCTATACTCCCTTTACATATTTATCATAACACAACGAGGCACACTTACACACCACATTAAAAACATAAAACCCTCACTCACACGAGAAAACATATTAATACTTATGCACCTCTTCCCCCTCCTCCTCCTAACCCTCAACCCTAACATCATTACTGGCTTTACTCCCTGTAAACATAGTTTAATCAAAACATTAGATTGTGAATCTAACAATAGAGGCTCG-AAACCTCTTGCTTACCGAGAAAGCC-CACAAGAACTGCTAACTCACTATCCCATGTATGACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAGCAATGTACACCACCATAGCCATTCTAACGCTAACCTCCCTAATTCCCCCCATTACAGCCACCCTTATTAACCCCAATAAAAAGAACTTATACCCGCACTACGTAAAAATGACCATTGCCTCTACCTTTATAATCAGCCTATTTCCCACAATAATATTCATGTGCACAGACCAAGAAACCATTATTTCAAACTGACACTGAACTGCAACCCAAACGCTAGAACTCTCCCTAAGCTT"));
        sequences.add(new Sequence("Macaca_fuscata","AAGCTTTTCCGGCGCAACCATCCTTATGATCGCTCACGGACTCACCTCTTCCATATATTTCTGCCTAGCCAATTCAAACTATGAACGCACTCACAACCGTACCATACTACTGTCCCGAGGACTTCAAATCCTACTTCCACTAACAGCCTTTTGATGATTAACAGCAAGCCTTACTAACCTTGCCCTACCCCCCACTATCAATCTACTAGGTGAACTCTTTGTAATCGCAACCTCATTCTCCTGATCCCATATCACCATTATGCTAACAGGACTTAACATATTAATTACGGCCCTCTACTCTCTCCACATATTCACTACAACACAACGAGGAACACTCACACATCACATAATCAACATAAAGCCCCCCTTCACACGAGAAAACACATTAATATTCATACACCTCGCTCCAATTATCCTTCTATCCCTCAACCCCAACATCATCCTGGGGTTTACCTCCTGTAGATATAGTTTAACTAAAACACTAGATTGTGAATCTAACCATAGAGACTCA-CCACCTCTTATTTACCGAGAAAACT-CGCAAGGACTGCTAACCCATGTACCCGTACCTAAAATTACGGTTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGACCTTAGGAGTCAAAAACATTGGTGCAACTCCAAATAAAAGTAATAATCATGCACACCCCCATCATTATAACAACCCTTATCTCCCTAACTCTCCCAATTTTTGCCACCCTCATCAACCCTTACAAAAAACGTCCATACCCAGATTACGTAAAAACAACCGTAATATATGCTTTCATCATCAGCCTCCCCTCAACAACTTTATTCATCTTCTCAAACCAAGAAACAACCATTTGGAGCTGACATTGAATAATGACCCAAACACTAGACCTAACGCTAAGCTT"));
        sequences.add(new Sequence("M_mulatta","AAGCTTTTCTGGCGCAACCATCCTCATGATTGCTCACGGACTCACCTCTTCCATATATTTCTGCCTAGCCAATTCAAACTATGAACGCACTCACAACCGTACCATACTACTGTCCCGGGGACTTCAAATCCTACTTCCACTAACAGCTTTCTGATGATTAACAGCAAGCCTTACTAACCTTGCCCTACCCCCCACTATCAACCTACTAGGTGAACTCTTTGTAATCGCGACCTCATTCTCCTGGTCCCATATCACCATTATATTAACAGGATTTAACATACTAATTACGGCCCTCTACTCCCTCCACATATTCACCACAACACAACGAGGAGCACTCACACATCACATAATCAACATAAAACCCCCCTTCACACGAGAAAACATATTAATATTCATACACCTCGCTCCAATCATCCTCCTATCTCTCAACCCCAACATCATCCTGGGGTTTACTTCCTGTAGATATAGTTTAACTAAAACATTAGATTGTGAATCTAACCATAGAGACTTA-CCACCTCTTATTTACCGAGAAAACT-CGCGAGGACTGCTAACCCATGTATCCGTACCTAAAATTACGGTTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGACCTTAGGAGTCAAAAATATTGGTGCAACTCCAAATAAAAGTAATAATCATGCACACCCCTATCATAATAACAACCCTTATCTCCCTAACTCTCCCAATTTTTGCCACCCTCATCAACCCTTACAAAAAACGTCCATACCCAGATTACGTAAAAACAACCGTAATATATGCTTTCATCATCAGCCTCCCCTCAACAACTTTATTCATCTTCTCAAACCAAGAAACAACCATTTGAAGCTGACATTGAATAATAACCCAAACACTAGACCTAACACTAAGCTT"));
        sequences.add(new Sequence("M_fascicularis","AAGCTTCTCCGGCGCAACCACCCTTATAATCGCCCACGGGCTCACCTCTTCCATGTATTTCTGCTTGGCCAATTCAAACTATGAGCGCACTCATAACCGTACCATACTACTATCCCGAGGACTTCAAATTCTACTTCCATTGACAGCCTTCTGATGACTCACAGCAAGCCTTACTAACCTTGCCCTACCCCCCACTATTAATCTACTAGGCGAACTCTTTGTAATCACAACTTCATTTTCCTGATCCCATATCACCATTGTGTTAACGGGCCTTAATATACTAATCACAGCCCTCTACTCTCTCCACATGTTCATTACAGTACAACGAGGAACACTCACACACCACATAATCAATATAAAACCCCCCTTCACACGAGAAAACATATTAATATTCATACACCTCGCTCCAATTATCCTTCTATCTCTCAACCCCAACATCATCCTGGGGTTTACCTCCTGTAAATATAGTTTAACTAAAACATTAGATTGTGAATCTAACTATAGAGGCCTA-CCACTTCTTATTTACCGAGAAAACT-CGCAAGGACTGCTAATCCATGCCTCCGTACTTAAAACTACGGTTTCCTCAACTTTTAAAGGATAACAGCTATCCATTGACCTTAGGAGTCAAAAACATTGGTGCAACTCCAAATAAAAGTAATAATCATGCACACCCCCATCATAATAACAACCCTCATCTCCCTGACCCTTCCAATTTTTGCCACCCTCACCAACCCCTATAAAAAACGTTCATACCCAGACTACGTAAAAACAACCGTAATATATGCTTTTATTACCAGTCTCCCCTCAACAACCCTATTCATCCTCTCAAACCAAGAAACAACCATTTGGAGTTGACATTGAATAACAACCCAAACATTAGACCTAACACTAAGCTT"));
        sequences.add(new Sequence("M_sylvanus","AAGCTTCTCCGGTGCAACTATCCTTATAGTTGCCCATGGACTCACCTCTTCCATATACTTCTGCTTGGCCAACTCAAACTACGAACGCACCCACAGCCGCATCATACTACTATCCCGAGGACTCCAAATCCTACTCCCACTAACAGCCTTCTGATGATTCACAGCAAGCCTTACTAATCTTGCTCTACCCTCCACTATTAATCTACTGGGCGAACTCTTCGTAATCGCAACCTCATTTTCCTGATCCCACATCACCATCATACTAACAGGACTGAACATACTAATTACAGCCCTCTACTCTCTTCACATATTCACCACAACACAACGAGGAGCGCTCACACACCACATAATTAACATAAAACCACCTTTCACACGAGAAAACATATTAATACTCATACACCTCGCTCCAATTATTCTTCTATCTCTTAACCCCAACATCATTCTAGGATTTACTTCCTGTAAATATAGTTTAATTAAAACATTAGACTGTGAATCTAACTATAGAAGCTTA-CCACTTCTTATTTACCGAGAAAACT-TGCAAGGACCGCTAATCCACACCTCCGTACTTAAAACTACGGTTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGCCTTAGGAGTCAAAAATATTGGTGCAACTCCAAATAAAAGTAATAATCATGTATACCCCCATCATAATAACAACTCTCATCTCCCTAACTCTTCCAATTTTCGCTACCCTTATCAACCCCAACAAAAAACACCTATATCCAAACTACGTAAAAACAGCCGTAATATATGCTTTCATTACCAGCCTCTCTTCAACAACTTTATATATATTCTTAAACCAAGAAACAATCATCTGAAGCTGGCACTGAATAATAACCCAAACACTAAGCCTAACATTAAGCTT"));
        sequences.add(new Sequence("Saimiri_sciureus","AAGCTTCACCGGCGCAATGATCCTAATAATCGCTCACGGGTTTACTTCGTCTATGCTATTCTGCCTAGCAAACTCAAATTACGAACGAATTCACAGCCGAACAATAACATTTACTCGAGGGCTCCAAACACTATTCCCGCTTATAGGCCTCTGATGACTCCTAGCAAATCTCGCTAACCTCGCCCTACCCACAGCTATTAATCTAGTAGGAGAATTACTCACAATCGTATCTTCCTTCTCTTGATCCAACTTTACTATTATATTCACAGGACTTAATATACTAATTACAGCACTCTACTCACTTCATATGTATGCCTCTACACAGCGAGGTCCACTTACATACAGCACCAGCAATATAAAACCAATATTTACACGAGAAAATACGCTAATATTTATACATATAACACCAATCCTCCTCCTTACCTTGAGCCCCAAGGTAATTATAGGACCCTCACCTTGTAATTATAGTTTAGCTAAAACATTAGATTGTGAATCTAATAATAGAAGAATA-TAACTTCTTAATTACCGAGAAAGTG-CGCAAGAACTGCTAATTCATGCTCCCAAGACTAACAACTTGGCTTCCTCAACTTTTAAAGGATAGTAGTTATCCATTGGTCTTAGGAGCCAAAAACATTGGTGCAACTCCAAATAAAAGTAATA---ATACACTTCTCCATCACTCTAATAACACTAATTAGCCTACTAGCGCCAATCCTAGCTACCCTCATTAACCCTAACAAAAGCACACTATACCCGTACTACGTAAAACTAGCCATCATCTACGCCCTCATTACCAGTACCTTATCTATAATATTCTTTATCCTTACAGGCCAAGAATCAATAATTTCAAACTGACACTGAATAACTATCCAAACCATCAAACTATCCCTAAGCTT"));
        Alignment alignment = new Alignment(sequences, 4, "nucleotide");
        
        // RecombinationGraph
        RecombinationGraph arg = new RecombinationGraph();
        ClusterTree tree = new ClusterTree();
        tree.initByName(
                "clusterType", "upgma",
                "taxa", alignment);
        
        arg.assignFrom(tree);
        arg.initAndValidate();
        
        // Site model:
        JukesCantor jc = new JukesCantor();
        jc.initByName();
        SiteModel siteModel = new SiteModel();
        siteModel.initByName(
                "substModel", jc);
        
        // Likelihood
        
        RecombinationGraphLikelihood argLikelihood = new RecombinationGraphLikelihood();
        argLikelihood.initByName(
                "data", alignment,
                "arg", arg,
                "siteModel", siteModel);
        
        arg.setEverythingDirty(true);
        System.out.println(argLikelihood.calculateLogP());
    }

}
