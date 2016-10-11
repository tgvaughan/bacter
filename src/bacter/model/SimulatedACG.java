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

import bacter.CFEventList;
import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.Locus;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Randomizer;
import feast.nexus.NexusBlock;
import feast.nexus.NexusBuilder;
import feast.nexus.TaxaBlock;
import feast.nexus.TreesBlock;

import java.io.FileNotFoundException;
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

    public Input<Double> rhoInput = new Input<>(
            "rho",
            "Conversion rate parameter.",
            Input.Validate.REQUIRED);
    
    public Input<Double> deltaInput = new Input<>(
            "delta",
            "Tract length parameter.",
            Input.Validate.REQUIRED);
    
    public Input<PopulationFunction> popFuncInput = new Input<>(
            "populationModel",
            "Demographic model to use.",
            Input.Validate.REQUIRED);

    public Input<Tree> clonalFrameInput = new Input<>(
            "clonalFrame",
            "Optional tree specifying fixed clonal frame.");

    public Input<String> outputFileNameInput = new Input<>(
            "outputFileName",
            "If provided, simulated ARG is additionally written to this file.");

    private double rho, delta;
    private PopulationFunction popFunc;

    public SimulatedACG() {
        m_taxonset.setRule(Input.Validate.REQUIRED);
    }
    
    @Override
    public void initAndValidate() {
        
        rho = rhoInput.get();
        delta = deltaInput.get();
        popFunc = popFuncInput.get();

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
            
            nexusBuilder.append(new TaxaBlock(m_taxonset.get()));
            
            nexusBuilder.append((new TreesBlock() {
                @Override
                public String getTreeString(Tree tree) {
                    return ((ConversionGraph)tree).getExtendedNewick();
                }
            }).addTree(this, "simulatedARG"));
            
            nexusBuilder.append(new NexusBlock() {

                @Override
                public String getBlockName() {
                    return "bacter";
                }

                @Override
                public List<String> getBlockLines() {
                    List<String> lines = new ArrayList<>();

                    String lociLine = "loci";
                    for (Locus locus : getLoci())
                        lociLine += " " + locus.getID() + ":" + locus.getSiteCount();
                    lines.add(lociLine);

                    lines.add("clonalframe_labeled " + root.toNewick());
                    lines.add("clonalframe_numbered " + root.toShortNewick(true));
                    for (Locus locus : getLoci()) {
                        for (Conversion conv : getConversions(locus)) {
                            lines.add("conversion node1=" + conv.getNode1().getNr()
                                    + " node2=" + conv.getNode2().getNr()
                                    + " site1=" + conv.getStartSite()
                                    + " site2=" + conv.getEndSite()
                                    + " locus=" + locus.getID());
                        }
                    }
                    
                    return lines;
                }
            });

            try (PrintStream pstream = new PrintStream(outputFileNameInput.get())) {
                nexusBuilder.write(pstream);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }
    }

    /**
     * Use coalescent model to simulate clonal frame.
     */
    private void simulateClonalFrame() {

        // Initialize leaf nodes
        List<Node> leafNodes = new ArrayList<>();
        for (int i=0; i<m_taxonset.get().getTaxonCount(); i++) {
            Node leaf = new Node();
            leaf.setNr(i);
            leaf.setID(m_taxonset.get().getTaxonId(i));
                        
            if (hasDateTrait())
                leaf.setHeight(getDateTrait().getValue(leaf.getID()));
            else
                leaf.setHeight(0.0);
            
            leafNodes.add(leaf);
        }
        
        // Create and sort list of inactive nodes
        List<Node> inactiveNodes = new ArrayList<>(leafNodes);
        Collections.sort(inactiveNodes, (Node n1, Node n2) -> {
            if (n1.getHeight()<n2.getHeight())
                return -1;
            
            if (n1.getHeight()>n2.getHeight())
                return 1;
            
            return 0;
        });
        
        List<Node> activeNodes = new ArrayList<>();
        
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
            (getTotalSequenceLength()+(delta-1.0)*getLoci().size()));

        // Generate conversions:
        for (int i=0; i<Nconv; i++) {
            // Choose alignment
            double u = Randomizer.nextDouble()*(getTotalSequenceLength()
                    + (delta-1.0)*getLoci().size());

            Locus affectedLocus = null;
            for (Locus locus : getLoci()) {
                if (u < locus.getSiteCount() + delta - 1.0) {
                    affectedLocus = locus;
                    break;
                } else
                    u -= locus.getSiteCount() + delta - 1.0;
            }

            if (affectedLocus == null)
                throw new IllegalStateException("Programmer error: " +
                        "locus choice loop fell through.");

            int startSite, endSite;
            if (u<delta) {
                startSite = 0;
            } else {
                startSite = (int)Math.ceil(u-delta);
            }
            endSite = startSite + (int)Randomizer.nextGeometric(1.0/delta);
            endSite = Math.min(endSite, affectedLocus.getSiteCount()-1);

            Conversion conv = new Conversion();
            conv.setLocus(affectedLocus);
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
