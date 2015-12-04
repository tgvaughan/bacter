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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Static methods for manipulating and interrogating lists of integers
 * representing integer ranges.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class IntRanges {

    /**
     * @param as1 range set argument 1
     * @param as2 range set argument 2
     * @return union between the two range sets
     */
    public static List<Integer> getUnion(List<Integer> as1, List<Integer> as2) {
        List<Integer> union = new ArrayList<>();
        int nextx, nexty;

        int i1 = 0, i2 = 0, ui = -2;
        int N1 = as1.size(), N2 = as2.size();
        while (i1 < N1 || i2 < N2) {
            if (i1 < N1) {
                if (i2 == N2 || as1.get(i1) < as2.get(i2)) {
                    nextx = as1.get(i1);
                    nexty = as1.get(i1 + 1);
                    i1 += 2;
                } else {
                    nextx = as2.get(i2);
                    nexty = as2.get(i2 + 1);
                    i2 += 2;
                }
            } else {
                nextx = as2.get(i2);
                nexty = as2.get(i2 + 1);
                i2 += 2;
            }

            if (ui < 0 || union.get(ui + 1) < nextx) {
                union.add(nextx);
                union.add(nexty);
                ui += 2;
            } else {
                if (union.get(ui + 1) < nexty)
                    union.set(ui + 1, nexty);
            }
        }

        return union;
    }

    /**
     * @param as1 range set argument 1
     * @param as2 range set argument 2
     * @return intersection between the two range sets
     */
    public static List<Integer> getIntersection(List<Integer> as1, List<Integer> as2) {
        List<Integer> intersection = new ArrayList<>();

        int i=0;
        int j=0;
        while (i<as1.size()) {

            while (j<as2.size()) {
                if (as2.get(j)>=as1.get(i+1))
                    break;

                if (as2.get(j+1)>as1.get(i)) {
                    if (as2.get(j)<as1.get(i))
                        intersection.add(as1.get(i));
                    else
                        intersection.add(as2.get(j));

                    if (as2.get(j+1)<=as1.get(i+1)) {
                        intersection.add(as2.get(j+1));
                    } else {
                        intersection.add(as1.get(i+1));
                        break;
                    }
                }

                j += 2;
            }
            i += 2;
        }

        return intersection;
    }

    /**
     * Partition ranges in as into ranges inside and outside of the
     * contiguous range [x,y].
     *
     * @param as range set
     * @param x left boundary of contiguous range
     * @param y right boundary of contiguous range
     * @param inside list to fill with inside ranges
     * @param outside list to fill with outside ranges
     */
    public static void partitionRanges(List<Integer> as, int x, int y,
                                       List<Integer> inside, List<Integer> outside) {

        int i=0;
        while (i<as.size() && as.get(i) < x)
            outside.add(as.get(i++));

        if (i%2==1) {
            outside.add(x);
            if (x<as.get(i))
                inside.add(x);
            else
                i += 1;
        }

        while (i<as.size() && as.get(i)<y)
            inside.add(as.get(i++));

        if (i%2==1) {
            inside.add(y);
            if (y<as.get(i))
                outside.add(y);
            else
                i += 1;
        }

        while (i<as.size())
            outside.add(as.get(i++));

    }

    /**
     * @param as range set
     * @return string representation of as.
     */
    public static String rangesToString(List<Integer> as) {
        String res = "";

        for (int i=0; i<as.size(); i+=2)
            res += " [" + as.get(i) + "," + as.get(i+1) + "]";

        return "{" + res + " }";
    }

    /**
     * Obtain the total number of sites included in a list of site ranges.
     *
     * @param as affected sites
     * @return total number of sites included
     */
    public static int getTotalSites(List<Integer> as) {
        int res = 0;

        for (int i=0; i<as.size(); i+=2) {
            res += as.get(i+1)-as.get(i);
        }

        return res;
    }

    /**
     * Test whether two range lists are disjoint.
     *
     * @param as1 list 1
     * @param as2 list 2
     * @return true iff all ranges are disjoint
     */
    public static boolean rangesAreDisjoint(List<Integer> as1, List<Integer> as2) {

        List<Integer> starts = new ArrayList<>();
        List<Integer> ends = new ArrayList<>();

        for (int i=0; i<as1.size(); i+=2) {
            starts.add(as1.get(i));
            ends.add(as1.get(i + 1));
        }

        for (int i=0; i<as2.size(); i+=2) {
            starts.add(as2.get(i));
            ends.add(as2.get(i+1));
        }

        Collections.sort(starts);
        Collections.sort(ends);

        for (int i=0; i<starts.size()-1; i++) {
            if (starts.get(i+1)<ends.get(i))
                return false;
        }

        return true;
    }

    /**
     * Clunky method which reads in a (sorted) range list from a string.
     * Assumes a string of the form "{ [1,3] [5,10] ... }".  Whitespace
     * is ignored and the curly braces are optional, but everything else
     * is assumed correct.  (No error checking is done!)
     *
     * @param string string to parse
     * @return resulting range list
     */
    public static List<Integer> fromString(String string) {
        List<Integer> result = new ArrayList<>();

        string = string.replaceAll("\\s","");

        if (string.startsWith("{"))
            string = string.substring(1, string.length() - 1);

        if (!string.isEmpty()) {
            string = string.substring(1, string.length() - 1);
            String[] pairs = string.split("\\]\\[");
            for (String pair : pairs) {
                String[] pairSplit = pair.split(",");
                result.add(Integer.parseInt(pairSplit[0]));
                result.add(Integer.parseInt(pairSplit[1]));
            }
        }

        return result;
    }

    /**
     * @param as1 first range
     * @param as2 second range
     * @return true if ranges as1 and as2 are equal.
     */
    public static boolean rangesEqual(List<Integer> as1, List<Integer> as2) {
        if (as1.size() != as2.size())
            return false;

        for (int i=0; i<as1.size(); i++) {
            if (!as1.get(i).equals(as2.get(i)))
                return false;
        }

        return true;
    }

    public static void main(String[] args) {
        List<Integer> as1 = new ArrayList<>();
        as1.add(3); as1.add(5);
        as1.add(9); as1.add(17);
        System.out.println("as1 = " + rangesToString(as1));

        List<Integer> as2 = new ArrayList<>();
        as2.add(1); as2.add(2);
        as2.add(5); as2.add(10);
        System.out.println("as2 = " + rangesToString(as2));
        List<Integer> as3 = getUnion(as1, as2);

        System.out.println("as3 = union(as1,as2) = " + rangesToString(as3));

        System.out.println("partitionRanges():");
        List<Integer> intersect = new ArrayList<>();
        List<Integer> difference = new ArrayList<>();
        int x = 4, y = 24;
        partitionRanges(as3, x, y, intersect, difference);
        System.out.format("as3 - [%d,%d] = %s\n", x, y, rangesToString(difference));
        System.out.format("as3 ^ [%d,%d] = %s\n", x, y, rangesToString(intersect));

        List<Integer> as4 = new ArrayList<>();
        as4.add(10); as4.add(20);
        as4.add(30); as4.add(40);
        List<Integer> as5 = new ArrayList<>();
        as5.add(1); as5.add(10);
        as5.add(20); as5.add(30);

        System.out.format("disjoint(%s, %s) = %s\n",
                rangesToString(as4),
                rangesToString(as5),
                rangesAreDisjoint(as4, as5));


        List<Integer> as6 = new ArrayList<>();
        as6.add(4); as6.add(24);

        System.out.format("getIntersection(as3, %s) = %s\n",
                rangesToString(as6),
                rangesToString(getIntersection(as3, as6)));
    }

}
