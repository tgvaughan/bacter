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
package bacter;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.StateNode;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import com.google.common.collect.Lists;
import feast.input.In;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Recombination graph based around the clonal frame.")
@Citation("Tim Vaughan, Alexei Drummod and Nigel French, "
        + "'Phylogenetic inference for bacteria in BEAST 2'. "
        + "In preparation.")
public class RecombinationGraph extends Tree {
    
    /**
     * Unlike Trees, Recombination graphs require an alignment (or at least
     * the length of an alignment) to be specified so that the regions of
     * the alignment affected by recombination events can be recorded.
     */
    public Input<Alignment> alignmentInput = In.create("alignment",
            "Sequence alignment corresponding to graph.");
    
    public Input<Integer> sequenceLengthInput = new In<Integer>(
            "sequenceLength",
            "Sequence length. Alternative to providing full alignment.")
            .setXOR(alignmentInput);
    
    public Input<String> fromStringInput = In.create("fromString",
            "Initialise ARG from string representation.");
    
    protected int sequenceLength;
    
    /**
     * List of recombinations on graph.
     */
    protected List<Recombination> recombs;
    protected List<Recombination> storedRecombs;
    
    /**
     * Class of events types on clonal frame.
     */
    public enum EventType {COALESCENCE, SAMPLE };

    /**
     * Class of events on clonal frame.
     */
    public class Event {
        EventType type;
        double t;
        int lineages;
        
        public Event(Node node) {
            if (node.isLeaf())
                type = EventType.SAMPLE;
            else
                type = EventType.COALESCENCE;
            
            t = node.getHeight();
        }
        
        public double getHeight() {
            return t;
        }
        
        public EventType getType() {
            return type;
        }
        
        /**
         * @return number of lineages _above_ this event.
         */
        public int getLineageCount() {
            return lineages;
        }
        
        @Override
        public String toString() {
            return "t: " + t + ", k: " + lineages + ", type: " + type;
        }
    }
    
    /**
     * List of events on clonal frame.
     */
    protected List<Event> cfEventList;
    protected boolean cfEventListDirty;
    
    
    @Override
    public void initAndValidate() throws Exception {

        recombs = Lists.newArrayList();
        storedRecombs = Lists.newArrayList();
        recombs.add(null); // Represents the clonal frame.
        
        if (alignmentInput.get() != null)
            sequenceLength = alignmentInput.get().getSiteCount();
        else
            sequenceLength = sequenceLengthInput.get();

        if (fromStringInput.get() != null) {
            fromString(fromStringInput.get());
        }
        
        cfEventList = Lists.newArrayList();
        cfEventListDirty = true;
        
        super.initAndValidate();
    }
    
    /**
     * Retrieve length of sequence, identifying bounds of recombination loci.
     * 
     * @return sequence length
     */
    public int getSequenceLength() {
        return sequenceLength;
    }
    
    /**
     * Add recombination to graph, ensuring recombination list
     * remains sorted.
     * 
     * @param recomb 
     */
    public void addRecombination(Recombination recomb) {
        startEditing();
        
        recomb.setRecombinationGraph(this);
        
        int i;
        for (i=1; i<recombs.size(); i++)
            if (recombs.get(i).startSite>recomb.startSite)
                break;
        
        recombs.add(i, recomb);
    }
    
    /**
     * Remove recombination from graph.
     * 
     * @param recomb recombination to remove.
     */
    public void deleteRecombination(Recombination recomb) {
        startEditing();
        
        if (recomb == null)
            throw new IllegalArgumentException("Cannot delete the clonal frame!");
        recombs.remove(recomb);
    }
    
    /**
     * Retrieve list of recombinations.
     * 
     * @return List of recombinations.
     */
    public List<Recombination> getRecombinations() {
        return recombs;
    }

    /**
     * Obtain number of recombination events.
     * 
     * @return Number of recombinations.
     */
    public int getNRecombs() {
        return recombs.size()-1;
    }
    
    /**
     * @return Total length of all edges in clonal frame.
     */
    public double getClonalFrameLength() {
        double length = 0.0;
        for (Node node : m_nodes) {
            if (node.isRoot())
                continue;
            length += node.getLength();
        }
        
        return length;
    }
    
    /**
     * Obtain the total number of sites corresponding to the clonal frame.
     * 
     * @return CF site count
     */
    public int getClonalFrameSiteCount() {
        int count = getSequenceLength();
        for (Recombination recomb : recombs) {
            if (recomb == null)
                continue;
            
            count -= recomb.getEndSite()-recomb.getStartSite() + 1;
        }
        return count;
    }

     /**
     * Check validity of recombinations.  Useful for probability densities
     * over the ARG to decide whether to return 0 based on an unphysical
     * state.
     * 
     * @return true if all recombinations are valid w.r.t. clonal frame.
     */
    public boolean isValid() {
        for (int ridx=1; ridx<recombs.size(); ridx++) {
            if (!recombs.get(ridx).isValid())
                return false;
            
            if (recombs.get(ridx).getStartSite()<0
                    || recombs.get(ridx).getStartSite()>=getSequenceLength()
                    || recombs.get(ridx).getEndSite()<0
                    || recombs.get(ridx).getEndSite()>=getSequenceLength())
                return false;
        }
        
        return true;
    }
    
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        
        for (Recombination recomb : getRecombinations()) {
            if (recomb == null)
                continue;
            sb.append(String.format("[&%d,%d,%s,%d,%d,%s] ",
                    recomb.node1.getNr(),
                    recomb.startSite,
                    String.valueOf(recomb.height1),
                    recomb.node2.getNr(),
                    recomb.endSite,
                    String.valueOf(recomb.height2)));
        }
        sb.append(super.toString());
        
        // Unfortunately, we must behave differently if we're being
        // called by toXML().
        StackTraceElement[] ste = Thread.currentThread().getStackTrace();
        if (ste[2].getMethodName().equals("toXML"))
            return sb.toString().replaceAll("&", "&amp;");
        else
            return sb.toString();
    }
    
    /**
     * Load ARG from string representation.  This is the same representation
     * used for XML state restoration.
     *
     * @param str
     */
    public void fromString(String str) {
        
        // Extract clonal frame and recombination components of string
        Pattern cfPattern = Pattern.compile("^[^\\(]*(\\(.*)$");
        Matcher cfMatcher = cfPattern.matcher(str);
        
        if (!cfMatcher.find())
            throw new RuntimeException("Error parsing ARG state string.");
        
        // Process clonal frame
        String sNewick = cfMatcher.group(cfMatcher.groupCount());
        try {
            TreeParser parser = new TreeParser();
            parser.thresholdInput.setValue(1e-10, parser);
            parser.offsetInput.setValue(0, parser);
            setRoot(parser.parseNewick(sNewick));
        } catch (Exception ex) {
            Logger.getLogger(RecombinationGraph.class.getName()).log(Level.SEVERE, null, ex);
        }

        initArrays();
        
        Pattern recombPattern = Pattern.compile("\\[&([^]]*)]");
        Matcher recombMatcher = recombPattern.matcher(str);
        
        // Process recombinations
        recombs.clear();
        recombs.add(null);
        
        while(recombMatcher.find()) {
            String [] elements = recombMatcher.group(1).split(",");
            
            Node node1 = getNode(Integer.parseInt(elements[0]));
            int startLocus = Integer.parseInt(elements[1]);
            double height1 = Double.parseDouble(elements[2]);
            
            Node node2 = getNode(Integer.parseInt(elements[3]));
            int endLocus = Integer.parseInt(elements[4]);
            double height2 = Double.parseDouble(elements[5]);

            Recombination recomb = new Recombination(
                    node1, height1,
                    node2, height2,
                    startLocus, endLocus);
            
            addRecombination(recomb);
        }
    }
    
    @Override
    public void fromXML(final org.w3c.dom.Node node) {
        fromString(node.getTextContent());
    }

    @Override
    public RecombinationGraph copy() {
        RecombinationGraph arg = new RecombinationGraph();
        
        arg.setID(getID());
        
        arg.index = index;
        arg.root = root.copy();
        arg.nodeCount = nodeCount;
        arg.internalNodeCount = internalNodeCount;
        arg.leafNodeCount = leafNodeCount;
        
        arg.recombs = Lists.newArrayList();
        arg.storedRecombs = Lists.newArrayList();
        
        for (Recombination recomb : recombs) {
            if (recomb == null)
                arg.recombs.add(null);
            else
                arg.recombs.add(recomb.getCopy());
        }
        arg.sequenceLength = sequenceLength;
        
        cfEventList = Lists.newArrayList();
        cfEventListDirty = true;

        return arg;
    }

    /**
     * Use another StateNode to configure this ARG.  If the other StateNode
     * is merely a tree, only the clonal frame is configured.
     * 
     * @param other StateNode used to configure ARG
     */
    @Override
    public void assignFrom(StateNode other) {
        super.assignFrom(other);
        
        if (other instanceof RecombinationGraph) {
            RecombinationGraph arg = (RecombinationGraph)other;
        
            recombs.clear();
            for (Recombination recomb : arg.recombs) {
                if (recomb == null)
                    recombs.add(null);
                else
                    recombs.add(recomb.getCopy());
            }
            sequenceLength = arg.sequenceLength;
        
            cfEventList = Lists.newArrayList();
            cfEventListDirty = true;
        }
    }
    
    @Override
    public void assignFromFragile(StateNode other) {
        super.assignFromFragile(other);
        
        RecombinationGraph arg = (RecombinationGraph)other;

        recombs.clear();
        storedRecombs.clear();
        for (Recombination recomb : arg.recombs) {
            if (recomb == null)
                recombs.add(null);
            else
                recombs.add(recomb.getCopy());
        }
        sequenceLength = arg.sequenceLength;
        
        cfEventList = Lists.newArrayList();
        cfEventListDirty = true;
    }
    
    /**
     * Obtain extended Newick representation of ARG.  If includeRegionInfo
     * is true, include Nexus metadata on hybrid leaf nodes describing the
     * alignment sites affected by the conversion event.
     * 
     * @param includeRegionInfo
     * @return Extended Newick string.
     */
    public String getExtendedNewick(boolean includeRegionInfo) {
        StringBuilder sb = new StringBuilder();
        
        sb.append(extendedNewickTraverse(root, includeRegionInfo));
        sb.append(";");
        
        return sb.toString();
    }
    
    private String extendedNewickTraverse(Node node, boolean includeRegionInfo) {
        StringBuilder sb = new StringBuilder();
        
        // Determine sequence of events along this node.
        class Event {
            boolean isArrival;
            double time;
            Recombination recomb;
            
            public Event(boolean isArrival, double time, Recombination recomb) {
                this.isArrival = isArrival;
                this.time = time;
                this.recomb = recomb;
            }
        }
        List<Event> events = new ArrayList<Event>();
        for (Recombination recomb : recombs) {
            if (recomb == null)
                continue;
            
            if (recomb.node1 == node)
                events.add(new Event(false, recomb.getHeight1(), recomb));
            if (recomb.node2 == node)
                events.add(new Event(true, recomb.getHeight2(), recomb));
        }
        
        // Sort events from oldest to youngest.
        Collections.sort(events, new Comparator<Event>() {
            @Override
            public int compare(Event e1, Event e2) {
                if (e1.time>e2.time)
                    return -1;
                else
                    return 1;
            }
        });

        // Process events.
        
        int cursor = 0;
        
        double lastTime;
        if (node.isRoot())
            lastTime = Double.POSITIVE_INFINITY;
        else
            lastTime = node.getParent().getHeight();

        for (Event event : events) {

            double thisLength;
            if (Double.isInfinite(lastTime))
                thisLength = 0.0;
            else
                thisLength = lastTime - event.time;
            
            if (event.isArrival) {
                String meta = !includeRegionInfo ? ""
                        : String.format("[&recomb=%d, region={%d,%d}]",
                                recombs.indexOf(event.recomb),
                                event.recomb.getStartSite(),
                                event.recomb.getEndSite());

                sb.insert(cursor, "(,#" + recombs.indexOf(event.recomb)
                        + meta
                        + ":" + (event.recomb.height2-event.recomb.height1)
                        + "):" + thisLength);
                cursor += 1;
            } else {
                sb.insert(cursor, "()#" + recombs.indexOf(event.recomb)
                        + ":" + thisLength);
                cursor += 1;
            }
            
            lastTime = event.time;
        }
        
        // Process this node and its children.

        if (!node.isLeaf()) {
            String subtree1 = extendedNewickTraverse(node.getChild(0),
                    includeRegionInfo);
            String subtree2 = extendedNewickTraverse(node.getChild(1),
                    includeRegionInfo);
            sb.insert(cursor, "(" + subtree1 + "," + subtree2 + ")");
            cursor += subtree1.length() + subtree2.length() + 3;
        }
        
        double thisLength;
        if (Double.isInfinite(lastTime))
            thisLength = 0.0;
        else
            thisLength = lastTime - node.getHeight();
        sb.insert(cursor, node.getNr() + ":" + thisLength);
        
        return sb.toString();
    }
    

    /**
     * Obtain ordered list of events that make up the clonal frame.  Used
     * for ARG probability density calculations and for various state proposal
     * operators.
     * 
     * @return List of events.
     */
    public List<Event> getCFEvents() {
        updateEvents();
        
        return cfEventList;
    }
    /**
     * Assemble sorted list of events on clonal frame and a map from nodes
     * to these events.
     */
    public void updateEvents() {
        if (!cfEventListDirty)
            return;
        
        cfEventList.clear();
        
        // Create event list
        for (Node node : getNodesAsArray()) {
            Event event = new Event(node);
            cfEventList.add(event);
        }
        
        // Sort events in increasing order of their times
        Collections.sort(cfEventList, new Comparator<Event>() {
            @Override
            public int compare(Event o1, Event o2) {
                if (o1.t<o2.t)
                    return -1;
                
                if (o2.t<o1.t)
                    return 1;
                
                return 0;
            }
        });
        
        // Compute lineage counts:
        int k=0;
        for (Event event : cfEventList) {
            if (event.type == EventType.SAMPLE)
                k += 1;
            else
                k -= 1;
            
            event.lineages = k;
        }
        
        cfEventListDirty = false;
    }
    
    /*
    * StateNode implementation
    */
    
    @Override
    protected void store () {
        super.store();
        
        storedRecombs.clear();

        for (Recombination recomb : recombs) {
            if (recomb == null)
                storedRecombs.add(null);
            else {
                Recombination recombCopy = new Recombination();
                
                recombCopy.setStartSite(recomb.getStartSite());
                recombCopy.setEndSite(recomb.getEndSite());
                recombCopy.setHeight1(recomb.getHeight1());
                recombCopy.setHeight2(recomb.getHeight2());
                
                recombCopy.setNode1(m_storedNodes[recomb.getNode1().getNr()]);
                recombCopy.setNode2(m_storedNodes[recomb.getNode2().getNr()]);

                recombCopy.setRecombinationGraph(this);
                
                storedRecombs.add(recombCopy);
            }
        }
    }
    
    @Override
    public void restore() {
        super.restore();
        
        List<Recombination> tmp = storedRecombs;
        storedRecombs = recombs;
        recombs = tmp;
        
        cfEventListDirty = true;
    }

    
    /**
     * Mark statenode as dirty if it belongs to a state.
     */
    public void startEditing() {
        if (state != null)
            startEditing(null);
        cfEventListDirty = true;
    }
    
    @Override
    public void startEditing(Operator operator) {
        super.startEditing(operator);
        cfEventListDirty = true;
    }
    
    /*
     * Loggable implementation.
     */
    
    @Override
    public void init(PrintStream out) throws Exception {

        out.println("#NEXUS");
        
        out.println("begin trees;");
        
        if (getTaxonset() != null) {
            StringBuilder translate = new StringBuilder("\ttranslate");
            for (int i=0; i<getLeafNodeCount(); i++) {
                if (i>0)
                    translate.append(",");
                translate.append(String.format(" %d %s",
                        i, getTaxonset().asStringList().get(i)));
            }
            out.println(translate + ";");
        }
        
    }
    
    
    @Override
    public void log(int nSample, PrintStream out) {
        RecombinationGraph arg = (RecombinationGraph) getCurrent();
        
        out.print(String.format("\ttree STATE_%d = [&R] %s",
                nSample, arg.getExtendedNewick(true)));
    }

    @Override
    public void close(PrintStream out) {
        out.println("end;");
    }
}
