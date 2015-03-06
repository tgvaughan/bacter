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
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import com.google.common.collect.Lists;
import feast.input.In;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Conversion graph based around the clonal frame.")
@Citation("Tim Vaughan, Alexei Drummod and Nigel French, "
        + "'Phylogenetic inference\n for bacteria in BEAST 2'. "
        + "In preparation.")
public class ConversionGraph extends Tree {
    
    /**
     * Unlike Trees, Conversion graphs require an alignment (or at least
 the length of an alignment) to be specified so that the regions of
 the alignment affected by recombination events can be recorded.
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
    protected List<Conversion> convs;
    protected List<Conversion> storedConvs;

    /**
     * Event and region lists.
     */
    protected RegionList regionList;
    protected CFEventList cfEventList;
    protected ACGEventList acgEventList;
    
    @Override
    public void initAndValidate() throws Exception {

        convs = Lists.newArrayList();
        storedConvs = Lists.newArrayList();
        
        if (alignmentInput.get() != null)
            sequenceLength = alignmentInput.get().getSiteCount();
        else
            sequenceLength = sequenceLengthInput.get();

        if (fromStringInput.get() != null) {
            fromString(fromStringInput.get());
        }

        regionList = new RegionList(this);
        cfEventList = new CFEventList(this);
        acgEventList = new ACGEventList(this);
        
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
     * @param conv
     */
    public void addConversion(Conversion conv) {
        startEditing(null);
        
        conv.setConversionGraph(this);
        
        int i;
        for (i=0; i<convs.size(); i++)
            if (convs.get(i).startSite>conv.startSite)
                break;
        
        convs.add(i, conv);
    }
    
    /**
     * Remove recombination from graph.
     * 
     * @param conv recombination to remove.
     */
    public void deleteConversion(Conversion conv) {
        startEditing(null);
        
        convs.remove(conv);
    }
    
    /**
     * Retrieve list of conversions.
     * 
     * @return List of conversions.
     */
    public List<Conversion> getConversions() {
        return convs;
    }

    /**
     * Obtain number of conversion events.
     * 
     * @return Number of conversions.
     */
    public int getConvCount() {
        return convs.size();
    }

    /**
     * Get list of contiguous regions having fixed marginal trees. 
     * 
     * @return list of regions
     */
    public List<Region> getRegions() {
        return regionList.getRegions();
    }

    /**
     * Obtain number of contiguous single-tree regions.
     * 
     * @return Number of regions.
     */
    public int getRegionCount() {
        return regionList.getRegions().size();
    }

    /**
     * Obtain ordered list of events that make up the clonal frame.  Used
     * for ACG probability density calculations and for various state proposal
     * operators.
     * 
     * @return List of events.
     */
    public List<CFEventList.Event> getCFEvents() {
        return cfEventList.getCFEvents();
    }

    /**
     * Obtain ordered list of events that make up the ACG.  Used
     * for ACG probability density calculations and for various state proposal
     * operators.
     * 
     * @return List of events.
     */
    public List<ACGEventList.Event> getACGEvents() {
        return acgEventList.getACGEvents();
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
        int count = 0;

        for (Region region : getRegions()) {
            if (region.isClonalFrame())
                count += region.leftBoundary - region.rightBoundary;
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
        for (Conversion conv : convs) {
            if (!conv.isValid()) {
                return false;
            }
            if (conv.getStartSite() < 0
                || conv.getStartSite() >= getSequenceLength()
                || conv.getEndSite() < 0
                || conv.getEndSite() >= getSequenceLength()) {
                return false;
            }
        }
        
        return true;
    }
    
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        
        for (Conversion conv : getConversions()) {
            sb.append(String.format("[&%d,%d,%s,%d,%d,%s] ",
                    conv.node1.getNr(),
                    conv.startSite,
                    String.valueOf(conv.height1),
                    conv.node2.getNr(),
                    conv.endSite,
                    String.valueOf(conv.height2)));
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
            Logger.getLogger(ConversionGraph.class.getName()).log(Level.SEVERE, null, ex);
        }

        initArrays();
        
        Pattern convPattern = Pattern.compile("\\[&([^]]*)]");
        Matcher convMatcher = convPattern.matcher(str);
        
        // Process recombinations
        convs.clear();
        
        while(convMatcher.find()) {
            String [] elements = convMatcher.group(1).split(",");
            
            Node node1 = getNode(Integer.parseInt(elements[0]));
            int startLocus = Integer.parseInt(elements[1]);
            double height1 = Double.parseDouble(elements[2]);
            
            Node node2 = getNode(Integer.parseInt(elements[3]));
            int endLocus = Integer.parseInt(elements[4]);
            double height2 = Double.parseDouble(elements[5]);

            Conversion recomb = new Conversion(
                    node1, height1,
                    node2, height2,
                    startLocus, endLocus);
            
            addConversion(recomb);
        }
    }
    
    @Override
    public void fromXML(final org.w3c.dom.Node node) {
        fromString(node.getTextContent());
    }

    @Override
    public ConversionGraph copy() {
        ConversionGraph arg = new ConversionGraph();
        
        arg.setID(getID());
        
        arg.index = index;
        arg.root = root.copy();
        arg.nodeCount = nodeCount;
        arg.internalNodeCount = internalNodeCount;
        arg.leafNodeCount = leafNodeCount;
        
        arg.convs = Lists.newArrayList();
        arg.storedConvs = Lists.newArrayList();
        
        for (Conversion recomb : convs) {
		arg.convs.add(recomb.getCopy());
        }
        arg.sequenceLength = sequenceLength;
        
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
        
        if (other instanceof ConversionGraph) {
            ConversionGraph arg = (ConversionGraph)other;
        
            convs.clear();
            for (Conversion recomb : arg.convs) {
                if (recomb == null)
                    convs.add(null);
                else
                    convs.add(recomb.getCopy());
            }
            sequenceLength = arg.sequenceLength;
        
            if (cfEventList == null)
                cfEventList = new CFEventList(this);
            cfEventList.makeDirty();

            if (acgEventList == null)
                acgEventList = new ACGEventList(this);
            acgEventList.makeDirty();

            if (regionList == null)
                regionList = new RegionList(this);
            regionList.makeDirty();
        }
    }
    
    @Override
    public void assignFromFragile(StateNode other) {
        super.assignFromFragile(other);
        
        ConversionGraph arg = (ConversionGraph)other;

        convs.clear();
        storedConvs.clear();
        for (Conversion recomb : arg.convs) {
            if (recomb == null)
                convs.add(null);
            else
                convs.add(recomb.getCopy());
        }
        sequenceLength = arg.sequenceLength;
        
        if (cfEventList == null)
            cfEventList = new CFEventList(null);
        cfEventList.makeDirty();

        if (acgEventList == null)
            acgEventList = new ACGEventList(this);
        acgEventList.makeDirty();
        
        if (regionList == null)
            regionList = new RegionList(this);
        regionList.makeDirty();
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
            Conversion recomb;
            
            public Event(boolean isArrival, double time, Conversion recomb) {
                this.isArrival = isArrival;
                this.time = time;
                this.recomb = recomb;
            }
        }
        List<Event> events = new ArrayList<>();
        for (Conversion recomb : convs) {
            if (recomb.node1 == node)
                events.add(new Event(false, recomb.getHeight1(), recomb));
            if (recomb.node2 == node)
                events.add(new Event(true, recomb.getHeight2(), recomb));
        }
        
        // Sort events from oldest to youngest.
        Collections.sort(events, (Event e1, Event e2) -> {
            if (e1.time>e2.time)
                return -1;
            else
                return 1;
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
                                convs.indexOf(event.recomb),
                                event.recomb.getStartSite(),
                                event.recomb.getEndSite());

                sb.insert(cursor, "(,#" + convs.indexOf(event.recomb)
                        + meta
                        + ":" + (event.recomb.height2-event.recomb.height1)
                        + "):" + thisLength);
                cursor += 1;
            } else {
                sb.insert(cursor, "()#" + convs.indexOf(event.recomb)
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


    
    /*
    * StateNode implementation
    */
    
    @Override
    protected void store () {
        super.store();
        
        storedConvs.clear();

        for (Conversion conv : convs) {
			Conversion convCopy = new Conversion();
	
			convCopy.setStartSite(conv.getStartSite());
			convCopy.setEndSite(conv.getEndSite());
			convCopy.setHeight1(conv.getHeight1());
			convCopy.setHeight2(conv.getHeight2());
			
			convCopy.setNode1(m_storedNodes[conv.getNode1().getNr()]);
			convCopy.setNode2(m_storedNodes[conv.getNode2().getNr()]);

			convCopy.setConversionGraph(this);
			
			storedConvs.add(convCopy);
        }
    }
    
    @Override
    public void restore() {
        super.restore();
        
        List<Conversion> tmp = storedConvs;
        storedConvs = convs;
        convs = tmp;

        cfEventList.makeDirty();
        acgEventList.makeDirty();
        regionList.makeDirty();
    }

    @Override
    public void startEditing(Operator operator) {
        if (state != null)
            super.startEditing(operator);

        if (cfEventList != null)
            cfEventList.makeDirty();

        if (acgEventList != null)
            acgEventList.makeDirty();

        if (regionList != null)
            regionList.makeDirty();
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
        ConversionGraph arg = (ConversionGraph) getCurrent();
        
        out.print(String.format("\ttree STATE_%d = [&R] %s",
                nSample, arg.getExtendedNewick(true)));
    }

    @Override
    public void close(PrintStream out) {
        out.println("end;");
    }
}
