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
import java.util.Objects;
import java.util.Set;

/**
 * Represents a contiguous region in which a single set of conversions is
 * active.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class Region {

    public int leftBoundary, rightBoundary;
    public Set<Conversion> activeConversions;

    public Region() {
        activeConversions = new HashSet<>();
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
}
