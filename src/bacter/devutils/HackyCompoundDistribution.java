package bacter.devutils;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.CompoundDistribution;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

public class HackyCompoundDistribution extends CompoundDistribution{
	
    public Input<TreeInterface> treeInput = new Input<TreeInterface>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

    public Input<SiteModelInterface.Base> siteModelInput = new Input<SiteModelInterface.Base>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);
    
    public Input<BranchRateModel.Base> branchRateModelInput = new Input<BranchRateModel.Base>("branchRateModel", "A model describing the rates on the branches of the beast.tree.", Validate.REQUIRED);
    
    public Input<Integer> stateCountInput = new Input<Integer>("stateCount", "A model describing the rates on the branches of the beast.tree.", Validate.REQUIRED);
	

	public Input<Buffer> bufferInput = new Input<Buffer>("buffer", "A model describing the rates on the branches of the beast.tree.", Validate.REQUIRED);
	
	TreeInterface tree;
	SiteModelInterface.Base siteModel; 
	BranchRateModel.Base branchRateModel;
	SubstitutionModel substitutionModel;
	
	double [][][] cfTransitionProbs;
	
	int stateCount;
	int categoryCount;
	
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
		tree = treeInput.get();
		branchRateModel = branchRateModelInput.get();
		siteModel = siteModelInput.get();
		substitutionModel = siteModel.getSubstitutionModel();
		
		stateCount = stateCountInput.get();
		categoryCount = siteModelInput.get().getCategoryCount();

		bufferInput.get().iniDimensions(tree.getNodeCount()-1, categoryCount, stateCount*stateCount);
		cfTransitionProbs = bufferInput.get().cfTransitionProbs;
		
	}
	
	@Override
    public double calculateLogP() {
		
		if(!ignoreInput.get()){
			preComputeTransitionProbs();
		}
		return super.calculateLogP();
    }

    /**
     * Pre-compute transition probabilities for CF edges.
     */
    void preComputeTransitionProbs() {
		
        for (int ni=0; ni<treeInput.get().getNodeCount()-1; ni++) {
            Node node = treeInput.get().getNode(ni);
            for (int ci = 0; ci < categoryCount; ci++) {
                double jointBranchRate = siteModelInput.get().getRateForCategory(ci, node) * branchRateModel.getRateForBranch(node);
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

}
