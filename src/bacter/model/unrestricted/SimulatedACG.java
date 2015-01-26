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

package bacter.model.unrestricted;

import bacter.CFEventList;
import bacter.Conversion;
import bacter.ConversionGraph;
import beast.core.Description;
import beast.core.Input;
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
import feast.nexus.TaxaBlock;
import feast.nexus.TreesBlock;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Simulates an ARG under the full ClonalOrigin model - can be used"
    + " for chain initialization or for sampler validation.")
public class SimulatedACG extends ConversionGraph {

    public Input<Double> rhoInput = new In<Double>("rho",
            "Conversion rate parameter.").setRequired();
    
    public Input<Double> deltaInput = new In<Double>("delta",
            "Tract length parameter.").setRequired();
    
    public Input<PopulationFunction> popFuncInput = new In<PopulationFunction>(
            "populationModel", "Demographic model to use.").setRequired();
    
    public Input<Integer> nTaxaInput = In.create("nTaxa",
            "Number of taxa to use in simulation. "
                    + "(Only use when alignment is unavailable.)");
    
    public Input<Tree> clonalFrameInput = In.create("clonalFrame",
            "Optional tree specifying fixed clonal frame.");

    public Input<String> outputFileNameInput = In.create("outputFileName",
            "If provided, simulated ARG is additionally written to this file.");

    private double rho, delta;
    private PopulationFunction popFunc;
    private int nTaxa;
    private TaxonSet taxonSet;
    
    public SimulatedACG() { };
    
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
                        List<Taxon> taxonList = new ArrayList<>();
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
        generateConversions();
        
        // Write output file
        if (outputFileNameInput.get() != null) {

            NexusBuilder nexusBuilder = new NexusBuilder();
            
            nexusBuilder.append(new TaxaBlock(taxonSet));
            
            nexusBuilder.append((new TreesBlock() {
                @Override
                public String getTreeString(Tree tree) {
                    return ((ConversionGraph)tree).getExtendedNewick(true);
                }
            }).addTree(this, "simulatedARG"));
            
            nexusBuilder.append(new NexusBlock() {

                @Override
                public String getBlockName() {
                    return "ARGBEAST";
                }

                @Override
                public List<String> getBlockLines() {
                    List<String> lines = new ArrayList<>();
                    lines.add("clonalframe " + root.toShortNewick(true));
                    for (int r = 0; r < getConvCount(); r++) {
                        Conversion recomb = getConversions().get(r);
                        lines.add("conversion node1=" + recomb.getNode1().getNr()
                        + " node2=" + recomb.getNode2().getNr()
                        + " site1=" + recomb.getStartSite()
                        + " site2=" + recomb.getEndSite());
                    }
                    
                    return lines;
                }
            });

            try (PrintStream pstream = new PrintStream(outputFileNameInput.get())) {
                nexusBuilder.write(pstream);
            }
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
        Collections.sort(inactiveNodes, (Node n1, Node n2) -> {
            if (n1.getHeight()<n2.getHeight())
                return -1;
            
            if (n2.getHeight()>n1.getHeight())
                return 1;
            
            return 0;
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
    
    private void generateConversions() {

        // Draw number of conversions:
        int Nconv = (int) Randomizer.nextPoisson(rho*getClonalFrameLength()*
            (getSequenceLength()+delta));

        // Generate conversions:
        for (int i=0; i<Nconv; i++) {
            int startSite, endSite;
            double u = Randomizer.nextDouble()*(delta + getSequenceLength());
            if (u<delta) {
                startSite = 0;
            } else {
                startSite = (int)(u-delta);
            }
            endSite = startSite + (int)Randomizer.nextGeometric(1.0/delta);
            endSite = Math.min(endSite, getSequenceLength()-1);

            Conversion conv = new Conversion();
            conv.setStartSite(startSite);
            conv.setEndSite(endSite);
            associateConversionWithCF(conv);
            addConversion(conv);
        }
    }
    
    /**
     * Associates recombination with the clonal frame, selecting points of
     * departure and coalescence.
     * 
     * @param conv recombination to associate
     */
    private void associateConversionWithCF(Conversion conv) {
    
        List<CFEventList.Event> eventList = getCFEvents();

        // Select departure point            
        double u = Randomizer.nextDouble()*getClonalFrameLength();
        
        boolean started = false;
        for (int eidx=0; eidx<eventList.size(); eidx++) {
            CFEventList.Event event = eventList.get(eidx);

            if (!started) {
                
                double interval = eventList.get(eidx+1).getHeight() - event.getHeight();
                
                if (u<interval*event.getLineageCount()) {
                    for (Node node : getNodesAsArray()) {
                        if (!node.isRoot()
                                && node.getHeight()<=event.getHeight()
                                && node.getParent().getHeight()>event.getHeight()) {
                            
                            if (u<interval) {
                                conv.setNode1(node);
                                conv.setHeight1(event.getHeight() + u);
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
                double t = Math.max(event.getHeight(), conv.getHeight1());

                double intervalArea;
                if (eidx<eventList.size()-1) {
                    intervalArea = popFunc.getIntegral(t, eventList.get(eidx+1).getHeight());
                } else
                    intervalArea = Double.POSITIVE_INFINITY;
                
                if (u<intervalArea*event.getLineageCount()) {
                    
                    // Fix height of attachment point

                    double tauEnd = popFunc.getIntensity(t) + u/event.getLineageCount();
                    double tEnd = popFunc.getInverseIntensity(tauEnd);
                    conv.setHeight2(tEnd);                    
                    
                    // Choose particular lineage to attach to
                    int nodeNumber = Randomizer.nextInt(event.getLineageCount());
                    for (Node node : getNodesAsArray()) {
                        if (node.getHeight()<=event.getHeight()
                                && (node.isRoot() || node.getParent().getHeight()>event.getHeight())) {
                            
                            if (nodeNumber == 0) {
                                conv.setNode2(node);
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
}
