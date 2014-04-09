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

package argbeast.operators;

import beast.core.Input;
import beast.core.Operator;
import feast.input.In;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ARGValidityWrapper extends RecombinationGraphOperator {
    
    public Input<Operator> operatorInput = new In<Operator>("operator",
            "Operator to follow validity check with.").setRequired();

    public ARGValidityWrapper() { }
    
    @Override
    public void initAndValidate() { }
    
    @Override
    public double proposal() {
        double logP = operatorInput.get().proposal();
        
        if (!argInput.get().isValid())
            return Double.NEGATIVE_INFINITY;
        
        return logP;
    }
    
}
