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

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import java.io.PrintStream;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class RecombinationGraphStats extends CalculationNode implements Loggable {

    public Input<RecombinationGraph> argInput = new Input <RecombinationGraph>(
            "arg", "Recombination graph to calculate summary statistics from.",
            Validate.REQUIRED);
    
    private RecombinationGraph arg;
    
    @Override
    public void initAndValidate() {
        arg = argInput.get();
    }
    
    private double getMeanTractLength() {
        
        if (arg.getNRecombs()<1)
            return Double.NaN;
        
        double mean = 0;
        for (Recombination recomb : arg.getRecombinations()) {
            if (recomb == null)
                continue;
            
            mean += recomb.getEndLocus()-recomb.getStartLocus()+1;
        }
        mean /= arg.getNRecombs();
        
        return mean;
    }
    
    private double getMeanInterTractLength() {
        
        if (arg.getNRecombs()<2)
            return Double.NaN;
        
        double mean = 0;
        for (int ridx=2; ridx<arg.getNRecombs(); ridx++) {
            mean += arg.getRecombinations().get(ridx).getStartLocus()
                    - arg.getRecombinations().get(ridx-1).getEndLocus() - 1;
        }
        mean /= arg.getNRecombs()-1;
        
        return mean;
    }
    
    private double getMeanEdgeLength() {
        
        if (arg.getNRecombs()<1)
            return Double.NaN;
        
        double mean = 0.0;
        for (Recombination recomb : arg.getRecombinations()) {
            if (recomb == null)
                continue;
            
            mean += recomb.getHeight2()-recomb.getHeight1();
        }
        mean /= arg.getNRecombs();
        
        return mean;
    }
    
    private double getMeanDepartureHeight() {
        if (arg.getNRecombs()<1)
            return Double.NaN;
        
        double mean = 0.0;
        for (Recombination recomb : arg.getRecombinations()) {
            if (recomb == null)
                continue;
            
            mean += recomb.getHeight1();
        }
        mean /= arg.getNRecombs();
        
        return mean;
    }
    
    @Override
    public void init(PrintStream out) throws Exception {
        String id = getID();
        if (id == null || id.matches("\\s*"))
            id = arg.getID();

        out.print(id + ".CFheight\t"
                + id + ".nRecomb\t"
                + id + ".meanTractLength\t"
                + id + ".meanInterTractLength\t"
                + id + ".meanEdgeLength\t"
                + id + ".meanDepartureHeight\t");
    }

    @Override
    public void log(int nSample, PrintStream out) {
        out.print(arg.getRoot().getHeight() + "\t"
                + arg.getNRecombs() + "\t"
                + getMeanTractLength() + "\t"
                + getMeanInterTractLength() + "\t"
                + getMeanEdgeLength() + "\t"
                + getMeanDepartureHeight() + "\t");
    }

    @Override
    public void close(PrintStream out) { }
    
}
