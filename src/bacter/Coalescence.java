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

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * Class of objects indicating lineages which are coalescing and at which
 * sites.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class Coalescence {
    public List<BitSet> descendantLeaves1, descendantLeaves2;
    public List<Integer> siteRanges;

    public Coalescence() {
        descendantLeaves1 = new ArrayList<>();
        descendantLeaves2 = new ArrayList<>();
        siteRanges = new ArrayList<>();
    }

    public Coalescence(String string) {

        siteRanges = new ArrayList<>();
        descendantLeaves1 = new ArrayList<>();
        descendantLeaves2 = new ArrayList<>();

        string = string.replaceAll("\\s+","");

        if (string.isEmpty())
            return;

        string = string.substring(1, string.length()-1);
        String[] split1 = string.replaceAll("\\s","").split("}\\[");

        for (String aSplit1 : split1) {
            String[] split2 = aSplit1.split("]\\{");
            String[] rangeStr = split2[0].split(",");

            siteRanges.add(Integer.parseInt(rangeStr[0]));
            siteRanges.add(Integer.parseInt(rangeStr[1]));

            String[] split3 = split2[1].split("\\}\\{");

            String[] bitStrA = split3[0].split(",");
            BitSet theseDecendentsA = new BitSet();
            for (String aBitStr : bitStrA)
                theseDecendentsA.set(Integer.parseInt(aBitStr));
            descendantLeaves1.add(theseDecendentsA);

            String[] bitStrB = split3[1].split(",");
            BitSet theseDecendentsB = new BitSet();
            for (String aBitStr : bitStrB)
                theseDecendentsB.set(Integer.parseInt(aBitStr));
            descendantLeaves2.add(theseDecendentsB);

        }
    }

    public int getIntervalCount() {
        return descendantLeaves1.size();
    }

    @Override
    public String toString() {
        String res = "";

        for (int i=0; i<getIntervalCount(); i++) {
            if (i>0)
                res += " ";

            res += "[" + siteRanges.get(2*i) + "," + siteRanges.get(2*i+1) + "]"
                    + descendantLeaves1.get(i).toString().replaceAll("\\s","")
                    + descendantLeaves2.get(i).toString().replaceAll("\\s","");
        }

        return res;
    }

    public void addInterval(int x, int y, BitSet dl1, BitSet dl2) {
        if (getIntervalCount()>0
                && x == siteRanges.get(siteRanges.size()-1)
                && ((descendantLeaves1.get(descendantLeaves1.size()-1).equals(dl1)
                && descendantLeaves2.get(descendantLeaves2.size()-1).equals(dl2))
                || (descendantLeaves1.get(descendantLeaves1.size()-1).equals(dl2)
                && descendantLeaves2.get(descendantLeaves2.size()-1).equals(dl1))))
            siteRanges.set(siteRanges.size()-1, y);
        else {
            siteRanges.add(x);
            siteRanges.add(y);
            descendantLeaves1.add(dl1);
            descendantLeaves2.add(dl2);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Coalescence that = (Coalescence) o;

        if (!siteRanges.equals(that.siteRanges))
            return false;

        for (int i=0; i<getIntervalCount(); i++) {
            if ((!descendantLeaves1.get(i).equals(that.descendantLeaves1.get(i))
                    && !descendantLeaves1.get(i).equals(that.descendantLeaves2.get(i)))
                    || (!descendantLeaves2.get(i).equals(that.descendantLeaves2.get(i))
                    && !descendantLeaves2.get(i).equals(that.descendantLeaves1.get(i))))
                return false;
        }

        return true;
    }

    @Override
    public int hashCode() {
        int result = siteRanges.hashCode();

        for (int i=0; i<getIntervalCount(); i++) {
            int dl1hash = descendantLeaves1.get(i).hashCode();
            int dl2hash = descendantLeaves2.get(i).hashCode();

            if (dl1hash < dl2hash)
                result = 31*(31 * result + dl1hash) + dl2hash;
            else
                result = 31*(31 * result + dl2hash) + dl1hash;
        }

        return result;
    }
}
