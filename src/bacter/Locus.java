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

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("A locus of a given length with which an ACG can be associated.")
public class Locus extends BEASTObject {

    public Input<Integer> siteCountInput = new Input<>(
            "siteCount",
            "Number of sites contained in alignments associated with this locus",
            Input.Validate.REQUIRED);

    protected int siteCount;

    public Locus() { }

    /**
     * Construct locus with given site count
     *
     * @param siteCount length of this locus
     */
    public Locus(int siteCount) {
        this.siteCount = siteCount;
    }

    @Override
    public void initAndValidate() throws Exception {
        siteCount = siteCountInput.get();
    }

    /**
     * @return number of sites associated with locus
     */
    public int getSiteCount() {
        return siteCount;
    }
}
