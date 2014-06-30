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

package argbeast.model;

import argbeast.Recombination;
import argbeast.RecombinationGraph;
import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.IntegerParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import feast.input.In;
import feast.nexus.NexusBlock;
import feast.nexus.NexusBuilder;
import feast.nexus.NexusWriter;
import feast.nexus.TaxaBlock;
import feast.nexus.TreesBlock;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Simulates an ARG - can be used for chain initialization or for "
        + "sampler validation.")
public class SimulatedRecombinationGraph extends RecombinationGraph implements StateNodeInitialiser {

    public Input<Double> rhoInput = new In<Double>("rho",
            "Recombination rate parameter.").setRequired();
    
    public Input<Double> deltaInput = new In<Double>("delta",
            "Tract length parameter.").setRequired();
    
    public Input<PopulationFunction> popFuncInput = new In<PopulationFunction>(
            "populationModel", "Demographic model to use.").setRequired();
    
    public Input<Integer> nTaxaInput = In.create("nTaxa",
            "Number of taxa to use in simulation. "
                    + "(Only use when alignment is unavailable.)");
    
    public Input<Tree> clonalFrameInput = In.create("clonalFrame",
            "Optional tree specifying fixed clonal frame.");
    
    public Input<IntegerParameter> mapInput = In.create(
            "recombinationMap", "Optional sequence of integers specifying "
                    + "sites affected by recombination events.  Fixes the "
                    + "total number of recombination events and the sites "
                    + "they affect, leaving only the clonal frame and "
                    + "recombinant edges to be simulated.");
    
    public Input<String> outputFileNameInput = In.create("outputFileName",
            "If provided, simulated ARG is additionally written to this file.");

    private double rho, delta;
    private PopulationFunction popFunc;
    private int nTaxa;
    private TaxonSet taxonSet;
    
    public SimulatedRecombinationGraph() { };
    
    @Override
    public void initAndValidate() throws Exception {
        
        rho = rhoInput.get();
        delta = deltaInput.get();
        popFunc = popFuncInput.get();
        
        if (alignmentInput.get() != null) {
            nTaxa = alignmentInput.get().getTaxonCount();
            taxonSet = new TaxonSet(alignmentInput.get());
        } else {
            if (clonalFrameInput.get() != null) {
                nTaxa = clonalFrameInput.get().getLeafNodeCount();
                taxonSet = clonalFrameInput.get().getTaxonset();
            } else {
                if (!m_traitList.get().isEmpty()) {
                    taxonSet = m_traitList.get().get(0).taxaInput.get();
                    nTaxa = taxonSet.getTaxonCount();
                } else {
                    if (nTaxaInput.get() != null) {
                        nTaxa = nTaxaInput.get();
                        List<Taxon> taxonList = new ArrayList<Taxon>();
                        for (int i=0; i<nTaxa; i++)
                            taxonList.add(new Taxon("t" + i));
                        taxonSet = new TaxonSet(taxonList);
                    } else
                        throw new IllegalArgumentException("Must specify nTaxa if"
                                + "neither alignment nor clonalFrame are used.");
                }
            }
        }
        m_taxonset.setValue(taxonSet, this);
        
        // Need to do this here as Tree.processTraits(), which is called
        // by hasDateTrait() and hence simulateClonalFrame(), expects a
        // tree with nodes.
        super.initAndValidate();
        
        if (clonalFrameInput.get() == null)
            simulateClonalFrame();
        else
            assignFromWithoutID(clonalFrameInput.get());
        
        // Need to do this here as this sets the tree object that the nodes
        // point to, so without it they point to the dummy tree created by
        // super.initAndValidate().
        initArrays();
        
        // Generate recombinations
        if (mapInput.get() == null)
            generateRecombinations();
        else {
            // Read recombination map directly from input
            if (mapInput.get().getDimension()%2 != 0)
                throw new IllegalArgumentException(
                        "Map must contain an even number of site indices.");
            
            for (int i=0; i<mapInput.get().getDimension()/2; i++) {
                int start = mapInput.get().getValue(i*2);
                int end = mapInput.get().getValue(i*2 + 1);
                if (end<start)
                    throw new IllegalArgumentException(
                            "Map site index pairs i,j must satisfy j>=i.");
                
                Recombination recomb = new Recombination();
                recomb.setStartSite(mapInput.get().getValue(i*2));
                recomb.setEndSite(mapInput.get().getValue(i*2 + 1));
                
                associateRecombinationWithCF(recomb);
                addRecombination(recomb);
            }
        }
        
        // Write output file
        if (outputFileNameInput.get() != null) {

            NexusBuilder nexusBuilder = new NexusBuilder();
            
            nexusBuilder.append(new TaxaBlock(taxonSet));
            
            nexusBuilder.append((new TreesBlock() {
                @Override
                public String getTreeString(Tree tree) {
                    return ((RecombinationGraph)tree).getExtendedNewick(true);
                }
            }).addTree(this, "simulatedARG"));
            
            nexusBuilder.append(new NexusBlock() {

                @Override
                public String getBlockName() {
                    return "ARGBEAST";
                }

                @Override
                public List<String> getBlockLines() {
                    List<String> lines = new ArrayList<String>();
                    lines.add("clonalframe " + root.toShortNewick(true));
                    for (int r = 1; r <= getNRecombs(); r++) {
                        Recombination recomb = getRecombinations().get(r);
                        lines.add("conversion node1=" + recomb.getNode1().getNr()
                        + " node2=" + recomb.getNode2().getNr()
                        + " site1=" + recomb.getStartSite()
                        + " site2=" + recomb.getEndSite());
                    }
                    
                    return lines;
                }
            });
            
            PrintStream pstream = new PrintStream(outputFileNameInput.get());
            nexusBuilder.write(pstream);
            pstream.close();
        }
    }
    
    /**
     * Use coalescent model to simulate clonal frame.
     */
    private void simulateClonalFrame() {

        // Initialize leaf nodes
        List<Node> leafNodes = Lists.newArrayList();
        for (int i=0; i<nTaxa; i++) {
            Node leaf = new Node();
            leaf.setNr(i);
            leaf.setID(taxonSet.asStringList().get(i));
                        
            if (hasDateTrait())
                leaf.setHeight(getDateTrait().getValue(i));
            else
                leaf.setHeight(0.0);
            
            leafNodes.add(leaf);
        }
        
        // Create and sort list of inactive nodes
        List<Node> inactiveNodes = Lists.newArrayList(leafNodes);
        Collections.sort(inactiveNodes, new Comparator<Node>() {
            @Override
            public int compare(Node n1, Node n2) {
                if (n1.getHeight()<n2.getHeight())
                    return -1;
                
                if (n2.getHeight()>n1.getHeight())
                    return 1;
                
                return 0;
            }
        });
        
        List<Node> activeNodes = Lists.newArrayList();
        
        double tau = 0.0;
        int nextNr = leafNodes.size();
        while (true) {
            
            // Calculate coalescence propensity
            int k = activeNodes.size();
            double chi = 0.5*k*(k-1);
            
            // Draw scaled coalescent time
            if (chi>0.0)
                tau += Randomizer.nextExponential(chi);
            else
                tau = Double.POSITIVE_INFINITY;
            
            // Convert to real time
            double t = popFunc.getInverseIntensity(tau);
            
            // If new time takes us past next sample time, insert that sample
            if (!inactiveNodes.isEmpty() && t>inactiveNodes.get(0).getHeight()) {
                Node nextActive = inactiveNodes.remove(0);
                activeNodes.add(nextActive);
                tau = popFunc.getIntensity(nextActive.getHeight());
                continue;
            }
            
            // Coalesce random pair of active nodes.
            Node node1 = activeNodes.remove(Randomizer.nextInt(k));
            Node node2 = activeNodes.remove(Randomizer.nextInt(k-1));
            
            Node parent = new Node();
            parent.addChild(node1);
            parent.addChild(node2);
            parent.setHeight(t);
            parent.setNr(nextNr++);
            
            activeNodes.add(parent);
            
            if (inactiveNodes.isEmpty() && activeNodes.size()<2)
                break;
        }
        
        // Remaining active node is root
        setRoot(activeNodes.get(0));
    }
    
    /**
     * Uses a discrete two-state Markov process to generate the list of
     * converted segments along the sequence.
     */
    private void generateRecombinations() {
        
        double pRec = 0.5*rho*getClonalFrameLength()/getSequenceLength();
        double pTractEnd = 1.0/delta;
        double p0cf = 1.0/(1.0 + pRec*delta);
        
        // Check for zero recombination rate (used sometimes for testing)
        if (pRec==0.0)
            return;
        
        int l; // next available convertible locus
        if (Randomizer.nextDouble()>p0cf) {
            Recombination recomb = new Recombination();
            recomb.setStartSite(0);
            int tractLength = (int)Randomizer.nextGeometric(pTractEnd);
            recomb.setEndSite(Math.min(tractLength, getSequenceLength()-1));
            associateRecombinationWithCF(recomb);
            addRecombination(recomb);
            
            l = tractLength + 1;
        } else
            l = 1;
        
        if (l>=getSequenceLength())
            return;
        
        while (true) {
            l += Randomizer.nextGeometric(pRec);

            if (l>=getSequenceLength())
                break;
            
            Recombination recomb = new Recombination();

            recomb.setStartSite(l);
            l += Randomizer.nextGeometric(pTractEnd);
            recomb.setEndSite(Math.min(l, getSequenceLength()-1));

            associateRecombinationWithCF(recomb);
            addRecombination(recomb);
            
            // The next site at which a conversion can begin
            l += 2;
            
            if (l>=getSequenceLength())
                break;
        }
    }
    
    /**
     * Associates recombination with the clonal frame, selecting points of
     * departure and coalescence.
     * 
     * @param recomb recombination to associate
     */
    private void associateRecombinationWithCF(Recombination recomb) {
    
        List<Event> eventList = getCFEvents();

        // Select departure point            
        double u = Randomizer.nextDouble()*getClonalFrameLength();
        
        boolean started = false;
        for (int eidx=0; eidx<eventList.size(); eidx++) {
            Event event = eventList.get(eidx);

            if (!started) {
                
                double interval = eventList.get(eidx+1).getHeight() - event.getHeight();
                
                if (u<interval*event.getLineageCount()) {
                    for (Node node : getNodesAsArray()) {
                        if (!node.isRoot()
                                && node.getHeight()<=event.getHeight()
                                && node.getParent().getHeight()>event.getHeight()) {
                            
                            if (u<interval) {
                                recomb.setNode1(node);
                                recomb.setHeight1(event.getHeight() + u);
                                break;
                            } else
                                u -= interval;
                            
                        }
                    }
                    started = true;
                    u = Randomizer.nextExponential(1.0);
                } else
                    u -= interval*event.getLineageCount();
            }
            
            if (started) {
                double t = Math.max(event.getHeight(), recomb.getHeight1());

                double intervalArea;
                if (eidx<eventList.size()-1) {
                    intervalArea = popFunc.getIntegral(t, eventList.get(eidx+1).getHeight());
                } else
                    intervalArea = Double.POSITIVE_INFINITY;
                
                if (u<intervalArea*event.getLineageCount()) {
                    
                    // Fix height of attachment point

                    double tauEnd = popFunc.getIntensity(t) + u/event.getLineageCount();
                    double tEnd = popFunc.getInverseIntensity(tauEnd);
                    recomb.setHeight2(tEnd);                    
                    
                    // Choose particular lineage to attach to
                    int nodeNumber = Randomizer.nextInt(event.getLineageCount());
                    for (Node node : getNodesAsArray()) {
                        if (node.getHeight()<=event.getHeight()
                                && (node.isRoot() || node.getParent().getHeight()>event.getHeight())) {
                            
                            if (nodeNumber == 0) {
                                recomb.setNode2(node);
                                break;
                            } else
                                nodeNumber -= 1;
                        }
                    }
                    break;
                } else
                    u -= intervalArea*event.getLineageCount();
            }
        }
    }
    
    @Override
    public void initStateNodes() throws Exception { }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(this);
    }
    
}
