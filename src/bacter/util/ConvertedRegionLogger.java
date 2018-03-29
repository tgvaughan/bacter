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

package bacter.util;

import bacter.Conversion;
import bacter.ConversionGraph;
import bacter.Locus;
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

    public Input<ConversionGraph> acgInput = new Input<>(
            "acg", "Conversion graph", Validate.REQUIRED);

    @Override
    public void initAndValidate() { }
    
    @Override
    public void init(PrintStream out) {
        final ConversionGraph arg = acgInput.get();

        String mainID = (getID() == null || getID().matches(("\\s*")))
                ? arg.getID() + ".converted"
                : getID();

        for (Locus locus : acgInput.get().getLoci())
                out.print(mainID + "." + locus.getID() + "\t");
    }

    @Override
    public void log(long nSample, PrintStream out) {

        for (Locus locus : acgInput.get().getLoci()) {
            if (acgInput.get().getConvCount(locus) == 0) {
                out.print("NA\t");
                return;
            }

            for (int r = 0; r < acgInput.get().getConvCount(locus); r++) {
                if (r > 0)
                    out.print(",");

                Conversion recomb = acgInput.get().getConversions(locus).get(r);
                out.print(recomb.getStartSite() + ":" + recomb.getEndSite());
            }
            out.print("\t");
        }
    }

    @Override
    public void close(PrintStream out) {
    }


}
