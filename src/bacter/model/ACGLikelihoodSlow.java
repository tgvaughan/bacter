/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
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
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("ACG likelihood which applies standard tree likelihood to" +
        " each marginal tree of a given ACG.  Very slow. For use in testing ONLY.")
public class ACGLikelihoodSlow extends GenericTreeLikelihood {

    public Input<Locus> locusInput = new Input<>(
            "locus",
            "Locus associated with alignment to evaluate probability of.",
            Input.Validate.REQUIRED);

    protected ConversionGraph acg;
    protected Locus locus;
    protected Alignment alignment;
    protected SubstitutionModel.Base substitutionModel;
    protected SiteModel.Base siteModel;

    @Override
    public void initAndValidate() throws Exception {
        if (treeInput.get() instanceof ConversionGraph)
            acg = (ConversionGraph) treeInput.get();
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

        siteModel = (SiteModel.Base) siteModelInput.get();
        substitutionModel = (SubstitutionModel.Base) siteModel.getSubstitutionModel();

        if (branchRateModelInput.get() != null &&
                !(branchRateModelInput.get() instanceof StrictClockModel))
            throw new IllegalArgumentException("ACGLikelihood currently only" +
                    "supports strict clock models.");

    }

    @Override
    public double calculateLogP() throws Exception {

        double logP = 0.0;
        for (Region region : acg.getRegions(locus)) {
            Alignment margAlign = createMarginalAlignment(alignment, acg, region);
            Tree margTree = new Tree(new MarginalTree(acg, region).getRoot());
            TreeLikelihood treeLikelihood = new TreeLikelihood();
            treeLikelihood.initByName(
                    "data", margAlign,
                    "tree", margTree,
                    "siteModel", siteModel);

            logP += treeLikelihood.calculateLogP();
        }

        return logP;

    }

    /**
     * Create Alignment object representing alignment of a region
     * corresponding to a single marginal tree.
     *
     * @param alignment
     * @param acg
     * @param region
     * @return
     * @throws Exception
     */
    public Alignment createMarginalAlignment(Alignment alignment,
                                             ConversionGraph acg, Region region) throws Exception {
        List<Sequence> sequences = new ArrayList<>();

        for (int leafIdx=0; leafIdx<alignment.getTaxonCount(); leafIdx++) {
            List<Integer> stateSequence;

            stateSequence = alignment.getCounts().get(leafIdx)
                    .subList(region.leftBoundary, region.rightBoundary);

            String taxonName = alignment.getTaxaNames().get(leafIdx);
            String charSequence = alignment.getDataType().state2string(stateSequence);

            sequences.add(new Sequence(taxonName, charSequence));
        }

        return new Alignment(sequences,
                alignment.getDataType().getTypeDescription());
    }
}
