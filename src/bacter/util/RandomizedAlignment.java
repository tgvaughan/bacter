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

import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class RandomizedAlignment extends Alignment {

    public Input<TaxonSet> taxonSetInput = new Input<>(
            "taxonSet",
            "Set of taxa to which alignment corresponds.",
            Input.Validate.REQUIRED);

    public Input<Integer> sequenceLengthInput = new Input<>(
            "siteCount",
            "Length of randomized alignment to generate.",
            Input.Validate.REQUIRED);

    public RandomizedAlignment() { }

    /**
     * Create a randomized alignment. Taxa are labeled 0 ... nTaxa-1.
     *
     * @param nTaxa Number of taxa in alignment.
     * @param length Number of sites in alignment.
     * @throws Exception
     */
    public RandomizedAlignment(int nTaxa, int length) throws Exception {

        List<Taxon> taxonList = new ArrayList<>();
        for (int i=0; i<nTaxa; i++) {
            taxonList.add(new Taxon(Integer.toString(i)));
        }

        taxonSetInput.setValue(new TaxonSet(taxonList), this);
        sequenceLengthInput.setValue(length, this);

        initAndValidate();
    }

    @Override
    public void initAndValidate() throws Exception {

        String[] alphabet = {"G", "T", "C", "A"};
        for (String taxonName : taxonSetInput.get().asStringList()) {
            StringBuilder sb = new StringBuilder();
            for (int i=0; i<sequenceLengthInput.get(); i++) {
                sb.append(alphabet[Randomizer.nextInt(4)]);
            }

            sequenceInput.setValue(new Sequence(taxonName, sb.toString()), this);
        }

        super.initAndValidate();
    }
}
