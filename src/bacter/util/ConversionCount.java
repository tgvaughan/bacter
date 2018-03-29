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
import beast.core.BEASTObject;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;

import java.io.PrintStream;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ConversionCount extends BEASTObject implements Loggable, Function {

    public Input<ConversionGraph> acgInput = new Input<>("acg",
            "Conversion graph", Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() {

    }

    // Loggable implementation

    @Override
    public void init(PrintStream out) {
        out.print(getID() + "\t");
    }

    @Override
    public void log(long nSample, PrintStream out) {
        out.print(acgInput.get().getTotalConvCount() + "\t");
    }

    @Override
    public void close(PrintStream out) {
    }

    // Function implementation

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue() {
        return acgInput.get().getTotalConvCount();
    }

    @Override
    public double getArrayValue(int iDim) {
        if (iDim != 0)
            throw new IndexOutOfBoundsException();

        return getArrayValue();
    }
}
