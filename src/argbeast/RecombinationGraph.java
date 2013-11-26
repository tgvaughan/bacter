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
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import com.google.common.collect.Lists;
import java.util.ArrayList;
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
public class RecombinationGraph extends Tree {
    
    /**
     * Unlike Trees, Recombination graphs require an alignment to be specified
     * so that the regions of the alignment affected by recombination events
     * can be recorded.
     */
    public Input<Alignment> alignmentInput = new Input<Alignment>("alignment",
            "Sequence alignment corresponding to graph.", Validate.REQUIRED);
    
    protected int sequenceLength;
    
    /**
     * List of recombinations on graph.
     */
    protected List<Recombination> recombs;
    
    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
        
        recombs = new ArrayList<Recombination>();
        recombs.add(null); // Represents the clonal frame.
        
        if (alignmentInput.get() != null)
            sequenceLength = alignmentInput.get().getSiteCount();
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
        int i;
        for (i=1; i<recombs.size(); i++)
            if (recombs.get(i).startLocus>recomb.startLocus)
                break;
        
        recombs.add(i, recomb);
    }
    
    /**
     * Remove recombination from graph.
     * 
     * @param recomb 
     */
    public void deleteRecombination(Recombination recomb) {
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
     * Retrieve root for marginal tree defined by recomb.
     * 
     * @param recomb
     * @return root node
     */
    public Node getMarginalRoot(Recombination recomb) {
        if (recomb==null)
            return root;
        
        if (recomb.node2.isRoot())
            return recomb.node1.getParent();
        else
            return getRoot();

    }
    
    /**
     * Retrieve parent of node in marginal tree defined by recomb.
     * 
     * @param node
     * @param recomb
     * @return node parent
     */
    public Node getMarginalParent(Node node, Recombination recomb) {
        if (recomb==null
                || recomb.node1==recomb.node2
                || recomb.node1.getParent()==recomb.node2)
            return node.getParent();
        
        if (node == recomb.node1 || node == recomb.node2)
            return recomb.node1.getParent();
        
        if (node == recomb.node1.getParent())
            return recomb.node2.getParent();
        
        if (node.getParent()==recomb.node1.getParent())
            return node.getParent().getParent();
        
        return node.getParent();
    }
    
    /**
     * Retrieve children of node in marginal tree defined by recomb.
     * 
     * @param node
     * @param recomb
     * @return node children
     */
    public List<Node> getMarginalChildren(Node node, Recombination recomb) {
        if (recomb==null
                || recomb.node1==recomb.node2
                || recomb.node1.getParent() == recomb.node2)
            return node.getChildren();
        
        List<Node> children = Lists.newArrayList();        
        
        // Recombination-shifted node
        if (node==recomb.node1.getParent()) {
            children.add(recomb.node1);
            children.add(recomb.node2);
            return children;
        }
        
        // Any other node:
        for (Node child : node.getChildren()) {
            if (child==recomb.node2)
                children.add(recomb.node1.getParent());
            else {
                if (child == recomb.node1.getParent()) {
                    if (child.getChild(0)==recomb.node1)
                        children.add(child.getChild(1));
                    else
                        children.add(child.getChild(0));
                } else
                    children.add(child);
            }
        }
        
        return children;
    }

    /**
     * Height of node in marginal tree defined by recomb.
     * 
     * @param node
     * @param recomb
     * @return node height
     */
    public double getMarginalNodeHeight(Node node, Recombination recomb) {
        if (recomb==null || recomb.node1==recomb.node2)
            return node.getHeight();

        if (node == recomb.node1.getParent())
            return recomb.height2;
        else
            return node.getHeight();
    }
    
    /**
     * Length of edge between node and parent in marginal tree defined
     * by recomb.
     * 
     * @param node
     * @param recomb
     * @return edge length (zero if node is root of marginal tree)
     */
    public double getMarginalBranchLength(Node node, Recombination recomb) {
        if (recomb==null || recomb.node1==recomb.node2)
            return node.getLength();
        
        Node parent = getMarginalParent(node, recomb);
        if (parent == null)
            return 0.0;
        else
            return getMarginalNodeHeight(parent, recomb)
                    -getMarginalNodeHeight(node, recomb);
    }
    
    /**
     * Determine whether node is root of marginal tree defined by recomb.
     * 
     * @param node
     * @param recomb
     * @return true if node is marginal root, false otherwise.
     */
    public boolean isNodeMarginalRoot(Node node, Recombination recomb) {
        if (recomb==null)
            return node.isRoot();
        
        return (getMarginalParent(node, recomb) == null);
    }

    /**
     * Determine whether node is leaf in marginal tree defined by recomb.
     * Just here for aesthetic completeness - all marginal trees share the
     * same leaves.
     * 
     * @param node
     * @param recomb
     * @return true if node is leaf, false otherwise.
     */
    public boolean isNodeMarginalLeaf(Node node, Recombination recomb) {
        return node.isLeaf();
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
            
            if (ridx>1) {
                if (recombs.get(ridx-1).getEndLocus()>recombs.get(ridx).getStartLocus()-2)
                    return false;
            } else {
                if (recombs.get(ridx).startLocus<0)
                    return false;
            }
            
            if (ridx==recombs.size()-1)
                if (recombs.get(ridx).getEndLocus()>getSequenceLength()-1)
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
                    recomb.startLocus,
                    String.valueOf(recomb.height1),
                    recomb.node2.getNr(),
                    recomb.endLocus,
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
    
    @Override
    public void fromXML(final org.w3c.dom.Node node) {
        String str = node.getTextContent();
        
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
    
    /**
     * Obtain extended Newick representation of ARG.
     * 
     * @return Extended Newick string.
     */
    public String getExtendedNewick() {
        StringBuilder sb = new StringBuilder();
        
        sb.append(extendedNewickTraverse(root));
        sb.append(";");
        
        return sb.toString();
    }
    
    private String extendedNewickTraverse(Node node) {
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
                sb.insert(cursor, "(,#" + recombs.indexOf(event.recomb)
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
            String subtree1 = extendedNewickTraverse(node.getChild(0));
            String subtree2 = extendedNewickTraverse(node.getChild(1));
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
     * Obtain Newick representation of marginal genealogy corresponding to
     * given recombination.
     * 
     * @param recomb Recombination object (null represents clonal frame)
     * @return Newick string
     */
    public String getMarginalNewick(Recombination recomb) {
        StringBuilder sb = new StringBuilder();
        
        sb.append(marginalNewickTraverse(getMarginalRoot(recomb), recomb));
        sb.append(";");

        return sb.toString();
    }
    
    /**
     * Generate Newick representation of marginal subtree below node
     * corresponding to given recombination.  Used only by getMarginalNewick().
     * 
     * @param node Node determining subtree
     * @param recomb Recombination resulting in marginal tree
     * @return Newick string
     */
    private String marginalNewickTraverse(Node node, Recombination recomb) {
        StringBuilder sb = new StringBuilder();
        
        if (!isNodeMarginalLeaf(node, recomb)) {
            sb.append("(");
            List<Node> children = getMarginalChildren(node, recomb);
            sb.append(marginalNewickTraverse(children.get(0), recomb));
            sb.append(",");
            sb.append(marginalNewickTraverse(children.get(1), recomb));
            sb.append(")");
        }
        
        sb.append(node.getNr());
        
        sb.append(":");
        if (isNodeMarginalRoot(node, recomb))
            sb.append("0.0");
        else
            sb.append(getMarginalBranchLength(node, recomb));
        
        return sb.toString();
    }
}
