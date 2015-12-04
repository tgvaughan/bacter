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

import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class IntRangesTest {

    @Test
    public void testUnion() throws Exception {
        List<Integer> as1 = IntRanges.fromString("[3,5] [9,17]");
        List<Integer> as2 = IntRanges.fromString("[1,2] [5,10]");

        List<Integer> union = IntRanges.getUnion(as1, as2);
        List<Integer> unionTruth = IntRanges.fromString("[1,2] [3,17]");

        assertTrue(IntRanges.rangesEqual(union, unionTruth));
    }

    @Test
    public void testPartition() throws Exception {
        List<Integer> as = IntRanges.fromString("[1,2] [3,17]");
        List<Integer> inside = new ArrayList<>();
        List<Integer> outside = new ArrayList<>();

        int x=4, y=24;

        IntRanges.partitionRanges(as, x, y, inside, outside);

        List<Integer> insideTruth = IntRanges.fromString("[4,17]");
        List<Integer> outsideTruth = IntRanges.fromString("[1,2] [3,4]");

        assertTrue(IntRanges.rangesEqual(inside, insideTruth));
        assertTrue(IntRanges.rangesEqual(outside, outsideTruth));
    }

    @Test
    public void testDisjoint() throws Exception {
        List<Integer> as1 = IntRanges.fromString("[3,5] [9,17]");
        List<Integer> as2 = IntRanges.fromString("[1,2] [5,10]");

    }
}
