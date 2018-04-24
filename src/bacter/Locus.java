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

package bacter;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("A locus of a given length with which an ACG can be associated.")
public class Locus extends BEASTObject {

    public Input<Integer> siteCountInput = new Input<>(
            "siteCount",
            "Number of sites contained in alignments associated with this locus");

    public Input<Alignment> alignmentInput = new Input<>(
            "alignment",
            "Initialize locus using this alignment.",
            Input.Validate.XOR, siteCountInput);

    public Input<Boolean> conversionsAllowedInput = new Input<>(
            "conversionsAllowed",
            "If false, no conversions will be allowed on this locus. (Default true.)",
            true);

    protected int siteCount;
    protected Alignment alignment;

    public Locus() { }

    /**
     * Construct locus corresponding to given alignment.
     *
     * @param alignment alignment object
     */
    public Locus(String name, Alignment alignment) {
        this.alignment = alignment;
        this.siteCount = alignment.getSiteCount();
        setID(name);
    }

    /**
     * Construct locus with given site count
     *
     * @param siteCount length of this locus
     */
    public Locus(String name, int siteCount) {
        this.alignment = null;
        this.siteCount = siteCount;
        setID(name);
    }

    @Override
    public void initAndValidate() {
        if (alignmentInput.get() != null) {
            alignment = alignmentInput.get();
            siteCount = alignment.getSiteCount();
        } else {
            alignment = null;
            siteCount = siteCountInput.get();
        }
    }

    /**
     * @return number of sites associated with locus
     */
    public int getSiteCount() {
        return siteCount;
    }

    /**
     * @return alignment associated with locus
     */
    public Alignment getAlignment() {
        return alignment;
    }

    /**
     * @return true if locus is associated with an alignment.
     */
    public boolean hasAlignment() {
        return alignment != null;
    }

    public boolean conversionsAllowed() {
        return conversionsAllowedInput.get();
    }
}
