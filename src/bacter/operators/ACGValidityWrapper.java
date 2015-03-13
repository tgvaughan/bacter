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

package bacter.operators;

import beast.core.Input;
import beast.core.Operator;
import feast.input.In;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ACGValidityWrapper extends ACGOperator {
    
    public Input<Operator> operatorInput = new In<Operator>("operator",
            "Operator to follow validity check with.").setRequired();

    public ACGValidityWrapper() {
        In.setOptional(m_pWeight);
    }

    @Override
    public void initAndValidate() throws Exception {
        m_pWeight.setValue(operatorInput.get().m_pWeight.get(), this);
        super.initAndValidate();
    }

    @Override
    public double proposal() {
        double logP = operatorInput.get().proposal();
        
        if (!acg.isValid())
            return Double.NEGATIVE_INFINITY;
        
        return logP;
    }

}