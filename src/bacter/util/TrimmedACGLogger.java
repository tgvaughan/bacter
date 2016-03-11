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

package bacter.util;

import bacter.ConversionGraph;
import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;

import java.io.PrintStream;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Logs ACGs with root-connecting conversions removed.")
public class TrimmedACGLogger extends CalculationNode implements Loggable {

    public Input<ConversionGraph> acgInput = new Input<>(
            "acg", "Conversion graph whose trimmed representation to log.",
            Validate.REQUIRED);

    @Override
    public void initAndValidate() { }
    
    @Override
    public void init(PrintStream out) {
       acgInput.get().init(out);
    }

    @Override
    public void log(int nSample, PrintStream out) {
        ConversionGraph arg = acgInput.get();

        out.print("tree STATE_" + nSample + " = ");
        out.print(arg.getTrimmedExtendedNewick());
    }

    @Override
    public void close(PrintStream out) {
        acgInput.get().close(out);
    }
    
}
