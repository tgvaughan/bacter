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
    public void testSplit1() {
        SiteAncestry a = new SiteAncestry("[200,400]{0} [600,800]{1}");

        SiteAncestry inside = new SiteAncestry();
        SiteAncestry outside = new SiteAncestry();
        a.split(250, 750, inside, outside);

        assertTrue(inside.equals(new SiteAncestry("[250,400]{0} [600,750]{1}")));
        assertTrue(outside.equals(new SiteAncestry("[200,250]{0} [750,800]{1}")));
    }


    @Test
    public void testSplit2() {
        SiteAncestry a = new SiteAncestry("[200,400]{0} [600,800]{1}");

        SiteAncestry inside = new SiteAncestry();
        SiteAncestry outside = new SiteAncestry();
        a.split(0, 100, inside, outside);

        assertTrue(inside.equals(new SiteAncestry("")));
        assertTrue(outside.equals(new SiteAncestry("[200,400]{0} [600,800]{1}")));
    }

    @Test
    public void testSplit3() {
        SiteAncestry a = new SiteAncestry("[200,400]{0} [600,800]{1}");

        SiteAncestry inside = new SiteAncestry();
        SiteAncestry outside = new SiteAncestry();
        a.split(500, 700, inside, outside);

        assertTrue(inside.equals(new SiteAncestry("[600, 700]{1}")));
        assertTrue(outside.equals(new SiteAncestry("[200,400]{0} [700,800]{1}")));
    }

    @Test
    public void testMerge1() {
        SiteAncestry a = new SiteAncestry("[0,100]{0} [200,300]{1}");
        SiteAncestry b = new SiteAncestry("[100,200]{0} [350,400]{1}");

        SiteAncestry union = new SiteAncestry();
        SiteAncestry coalescence = new SiteAncestry();
        a.merge(b, coalescence, union);

        assertTrue(union.equals(new SiteAncestry("[0,200]{0} [200,300]{1} [350,400]{1}")));
    }

    @Test
    public void testMerge2() {
        SiteAncestry a = new SiteAncestry("[0,100]{0} [120,150]{0} [200,300]{1}");
        SiteAncestry b = new SiteAncestry("[100,200]{1} [250,400]{2}");

        SiteAncestry union = new SiteAncestry();
        SiteAncestry coalescence = new SiteAncestry();
        a.merge(b, coalescence, union);

        assertTrue(union.equals(new SiteAncestry("[0,100]{0} [100,120]{1} " +
                "[120,150]{0,1} [150,250]{1} [250,300]{1,2} [300,400]{2}")));
    }
}
