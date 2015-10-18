package bacter;

import com.google.common.collect.Range;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class AffectedSiteList {

    ConversionGraph acg;
    Map<Conversion, List<Integer>> affectedSites;

    public AffectedSiteList(ConversionGraph acg) {
        this.acg = acg;
    }

    public void computeAffectedSites() {

    }

    static List<Integer> getUnion(List<Integer> as1, List<Integer> as2) {
        List<Integer> union = new ArrayList<>();
        int nextx, nexty;

        int i1 = 0, i2 = 0, ui = -2;
        int N1 = as1.size(), N2 = as2.size();
        while (i1 < N1 || i2 < N2) {
            if (i1<N1) {
                if (i2==N2 || as1.get(i1) < as2.get(i2)) {
                    nextx = as1.get(i1);
                    nexty = as1.get(i1+1);
                    i1 += 2;
                } else {
                    nextx = as2.get(i2);
                    nexty = as2.get(i2+1);
                    i2 += 2;
                }
            } else {
                nextx = as2.get(i2);
                nexty = as2.get(i2+1);
                i2 += 2;
            }

            if (ui<0 || union.get(ui+1)<nextx) {
                union.add(nextx);
                union.add(nexty);
                ui += 2;
            } else {
                if (union.get(ui+1)<nexty)
                    union.set(ui+1, nexty);
            }
        }

        return union;
    }

    public static List<Integer> getIntersection(List<Integer> as, int x, int y) {
        List<Integer> intersect = new ArrayList<>();

        // Early exit for empty range
        if (y==x)
            return intersect;

        int ix = Collections.binarySearch(as, x);
        if (ix<0)
            ix = -(ix+1)-1;

        int iy = Collections.binarySearch(as, y);
        if (iy<0)
            iy = -(iy+1)-1;


        if (ix == iy && ix%2 == 0) {
            intersect.add(x);
            intersect.add(y);

            return intersect;
        }

        for (int i=ix; i<=iy; i++) {
            if (i%2 == 0) {
                if (as.get(i)<x) {
                    intersect.add(x);
                } else {
                    intersect.add(as.get(i));
                }
            } else {
                intersect.add(as.get(i));
            }
        }
        if (ix%2 == 0) {
            intersect.add(x);
            ix += 1;
        }

        if (iy>ix)
            intersect.addAll(as.subList(ix, iy));

        if (iy%2 == 1) {
            intersect.add(as.get(iy));
        } else if (y>as.get(iy)) {
            intersect.add(as.get(iy));
            intersect.add(y);
        }


        return intersect;
    }

    static String rangesToString(List<Integer> bounds) {
        String res = "";

        for (int i=0; i<bounds.size(); i+=2)
            res += " [" + bounds.get(i) + "," + bounds.get(i+1) + "]";

        return "{" + res + " }";
    }

    public static void main(String[] args) {
        List<Integer> as1 = new ArrayList<>();
        as1.add(3); as1.add(5);
        as1.add(9); as1.add(17);
        System.out.println("as1 = " + rangesToString(as1));

        List<Integer> as2 = new ArrayList<>();
        as2.add(23); as2.add(25);
        System.out.println("as2 = " + rangesToString(as2));
        List<Integer> as3 = getUnion(as1, as2);

        System.out.println("as3 = union(as1,as2) = " + rangesToString(as3));

        List<Integer> as4 = getIntersection(as3, 8, 19);
        System.out.println("as4 = " + rangesToString(as4));


    }

}
