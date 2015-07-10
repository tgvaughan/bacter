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
package bacter;

import java.util.HashSet;
import java.util.Set;

/**
 * Represents a contiguous region in which a single set of conversions is
 * active.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class Region {

    final public int leftBoundary, rightBoundary;
    final public Set<Conversion> activeConversions;

    /**
     * Construct new region.
     *
     * @param leftBoundary leftmost boundary
     * @param rightBoundary rightmost boundary
     * @param activeConversions set of active conversions
     */
    public Region(int leftBoundary, int rightBoundary, Set<Conversion> activeConversions) {
        this.leftBoundary = leftBoundary;
        this.rightBoundary = rightBoundary;
        this.activeConversions = new HashSet<>(activeConversions);
    }

    /**
     * Construct new clonal frame region.
     *
     * @param leftBoundary leftmost boundary
     * @param rightBoundary rightmost boundary
     */
    public Region(int leftBoundary, int rightBoundary) {
        this.leftBoundary = leftBoundary;
        this.rightBoundary = rightBoundary;
        this.activeConversions = new HashSet<>();
    }


    public int getRegionLength() {
        return rightBoundary - leftBoundary;
    }

    public boolean isClonalFrame() {
        return activeConversions.isEmpty();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(leftBoundary).append("-").append(rightBoundary-1);
        if (isClonalFrame())
            sb.append(" (CF)");
        else
            sb.append(" (")
                .append(activeConversions.size())
                .append(" active conversions)");

        return sb.toString();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Region region = (Region) o;

        if (leftBoundary != region.leftBoundary) return false;
        if (rightBoundary != region.rightBoundary) return false;
        return activeConversions.equals(region.activeConversions);

    }

    @Override
    public int hashCode() {
        int result = leftBoundary;
        result = 31 * result + rightBoundary;
        result = 31 * result + activeConversions.hashCode();
        return result;
    }
}
