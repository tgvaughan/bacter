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

package argbeast.util;

import argbeast.Recombination;
import argbeast.RecombinationGraph;
import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import java.io.PrintStream;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("")
public class ConvertedRegionLogger extends BEASTObject implements Loggable {

    public Input<RecombinationGraph> argInput = new Input<RecombinationGraph>(
            "arg", "Recombination graph", Validate.REQUIRED);

    @Override
    public void initAndValidate() { }
    
    @Override
    public void init(PrintStream out) throws Exception {
        final RecombinationGraph arg = argInput.get();
        if (getID() == null || getID().matches("\\s*")) {
            out.print(arg.getID() + ".converted\t");
        } else {
            out.print(getID() + "\t");
        }
    }

    @Override
    public void log(int nSample, PrintStream out) {
        
        if (argInput.get().getNRecombs()==0) {
            out.print("NA\t");
            return;
        }
        
        for (int r=1; r<=argInput.get().getNRecombs(); r++) {
            if (r>1)
                out.print(",");

            Recombination recomb = argInput.get().getRecombinations().get(r);
            out.print(recomb.getStartSite() + ":" + recomb.getEndSite());
        }
        out.print("\t");
    }

    @Override
    public void close(PrintStream out) {
    }


}
