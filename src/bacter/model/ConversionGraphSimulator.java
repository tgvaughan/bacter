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

import beast.core.Description;
import beast.core.Input;
import beast.core.Logger;
import feast.input.In;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Simulates an ARG - can be used for chain initialization or for "
        + "sampler validation.")
public class ConversionGraphSimulator extends beast.core.Runnable {

    public Input<SimulatedConversionGraph> simARGInput =
            new In<SimulatedConversionGraph>("simARG",
                    "Simulated recombination graph.").setRequired();

    public Input<Integer> nSimsInput = new In<Integer>("nSims",
            "Number of ARGs to simulate.").setRequired();
    
    public Input<List<Logger>> loggersInput = new In<List<Logger>>("logger",
            "Logger used to write results to screen or disk.")
            .setDefault(new ArrayList<>());    
    
    @Override
    public void initAndValidate() { }

    @Override
    public void run() throws Exception {

        // Initialise loggers
        for (Logger logger : loggersInput.get()) {
            logger.init();
        }

        for (int i=0; i<nSimsInput.get(); i++) {
            
            if (i>0)
                simARGInput.get().initAndValidate();
            
            // Log state
            for (Logger logger : loggersInput.get()) {
                logger.log(i);
            }
        }
        
        // Finalize loggers
        for (Logger logger : loggersInput.get()) {
            logger.close();
        }
    }
}
