/*
 * Copyright (C) 2016 Tim Vaughan <tgvaughan@gmail.com>
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
import beast.core.Description;
import beast.core.Input;
import beast.evolution.sitemodel.SiteModel;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Simulates an alignment down an existing ACG.")
public class AlignmentSimulator extends beast.core.Runnable {

    public Input<ConversionGraph> acgInput = new Input<>(
            "acg",
            "Conversions graph down which to simulate evolution.",
            Input.Validate.REQUIRED);

    public Input<SiteModel> siteModelInput = new Input<>(
            "siteModel",
            "site model for leafs in the beast.tree",
            Input.Validate.REQUIRED);

    public Input<String> outputFileNameInput = new Input<>(
            "outputFileName",
            "If provided, simulated alignment is additionally written to this file.");

    public Input<Boolean> useNexusInput = new Input<>(
            "useNexus",
            "Use Nexus instead of FASTA format to write alignment file.",
            false);

    public Input<Locus> locusInput = new Input<>("locus",
            "Locus for which alignment will be simulated.");

    @Override
    public void initAndValidate() { }

    @Override
    public void run() throws Exception {
        SimulatedAlignment alignment = new SimulatedAlignment();
        alignment.initByName(
                "acg", acgInput.get(),
                "siteModel", siteModelInput.get(),
                "outputFileName", outputFileNameInput.get(),
                "useNexus", useNexusInput.get(),
                "locus", locusInput.get());
    }
}
