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

package bacter.devutils;

import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.MCMC;
import beast.base.inference.State;
import beast.base.parser.XMLParser;

import java.io.File;

/**
 * Computes density of given state under given model
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class DensityCalculator {

    public static void main(String[] args) {

        if (args.length < 2) {
            System.out.println("Usage: DensityCalculator BEAST_xml_file state_file");
            System.exit(0);
        }

        XMLParser parser = new XMLParser();
        beast.base.inference.Runnable runnable = null;
        try {
            runnable = parser.parseFile(new File(args[0]));
        } catch (Exception e) {
            System.out.println("Encountered error while loading/parsing XML file.");
            e.printStackTrace();
        }

        if (runnable == null)
            System.exit(1);

        if (!(runnable instanceof MCMC)) {
            System.out.println("XML file does not seem to describe an MCMC analysis.");
            System.exit(1);
        }

        MCMC mcmc = (MCMC)runnable;
        Distribution posterior = mcmc.posteriorInput.get();

        State state = mcmc.startStateInput.get();

        state.setStateFileName(args[1]);
        try {
            state.restoreFromFile();
        } catch (Exception e) {
            System.out.println("Error reading state from file.");
            e.printStackTrace();
            System.exit(1);
        }

        try {
            state.robustlyCalcPosterior(posterior);
        } catch (Exception e) {
            e.printStackTrace();
        }

        reportTargetDensities(posterior);
    }

    protected static void reportTargetDensities(final Distribution distr) {
        System.out.println(distr.getID() + ": " + distr.getCurrentLogP());
        if (distr instanceof CompoundDistribution) {
            for (final Distribution childDistr : ((CompoundDistribution) distr).pDistributions.get()) {
                reportTargetDensities(childDistr);
            }
        }
    }
}
