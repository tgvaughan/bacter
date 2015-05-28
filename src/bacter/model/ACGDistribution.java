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

package bacter.model;

import bacter.ConversionGraph;
import beast.core.Distribution;
import beast.core.Input;
import feast.input.In;

/**
 * Abstract class of distributions over ConversionGraphs. Concrete
 * implementations should override calculateACGLogP rather than
 * calculateLogP.  This allows optional validity testing to be performed.
 * 
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public abstract class ACGDistribution extends Distribution {

    public Input<ConversionGraph> acgInput = new In<ConversionGraph>(
            "acg", "Recombination graph.").setRequired();

    public Input<Boolean> checkValidityInput = new In<Boolean>("checkValidity",
            "Explicitly check validity of conversion graph.").setDefault(false);

    protected ConversionGraph acg;

    @Override
    public void initAndValidate() throws Exception {
        acg = acgInput.get();
    }

    @Override
    public final double calculateLogP() throws Exception {
        if (checkValidityInput.get() && acg.isInvalid()) {
            logP = Double.NEGATIVE_INFINITY;
        } else {
            logP = calculateACGLogP();
        }

        return logP;
    }

    public abstract double calculateACGLogP();
    
}
