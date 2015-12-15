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

import org.junit.Test;

import static org.junit.Assert.assertTrue;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SiteAncestryTest {

    @Test
    public void testSplit() {
        SiteAncestry a = new SiteAncestry("[0,400]{0} [600,1000]{1}");

        SiteAncestry inside = new SiteAncestry();
        SiteAncestry outside = new SiteAncestry();
        a.split(250, 750, inside, outside);

        assertTrue(inside.equals(new SiteAncestry("[250,400]{0} [600,750]{1}")));
        assertTrue(outside.equals(new SiteAncestry("[0,250]{0} [750,1000]{1}")));
    }
}
