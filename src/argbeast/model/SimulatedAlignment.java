/*
 * Copyright (C) 2014 Tim Vaughan <tgvaughan@gmail.com>
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

import argbeast.RecombinationGraph;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.sitemodel.SiteModel;
import feast.input.In;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("An alignment produced by simulating sequence evolution down an ARG.")
public class SimulatedAlignment extends Alignment {
    
    public Input<RecombinationGraph> argInput = new In<RecombinationGraph>(
            "ARG",
            "Recombination graph down which to simulate evolution.")
            .setRequired();
    
    public Input<SiteModel.Base> siteModelInput = new In<SiteModel.Base>(
            "siteModel",
            "site model for leafs in the beast.tree")
            .setRequired();
    
    public Input<BranchRateModel.Base> branchRateModelInput = In.create(
            "branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");
    
    public Input<String> outputFileNameInput = In.create(
            "outputFileName",
            "If provided, simulated alignment is additionally written to this file.");    
    
    public SimulatedAlignment() {
        
        // Override the sequence input requirement.
        In.setOptional(sequenceInput);
    }

    @Override
    public void initAndValidate() throws Exception {

        super.initAndValidate();
    }
    
}
