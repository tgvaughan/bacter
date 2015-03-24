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

package bacter.operators.unrestricted;

import bacter.Conversion;
import beast.core.Description;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

@Description("This operator replaces a single conversion with a pair of" +
        " conversions, causing the lineage to temporarily visit another" +
        " clonal frame edge.")
public class AddRemoveDetour extends ConversionCreationOperator {

        @Override
        public double proposal() {
                if (Randomizer.nextBoolean())
                        return addDetour();
                else
                        return removeDetour();
        }

        /**
         * Detour addition variant of move.
         * @return log HGF
         */
        private double addDetour() {
                double logHGF = 0.0;

                // Select conversion at random
                Conversion conv = acg.getConversions().get(
                        Randomizer.nextInt(acg.getConvCount()));
                logHGF -= Math.log(1.0/acg.getConvCount());

                // Select detour times:
                double t1 = conv.getHeight1()
                        + Randomizer.nextDouble()*(conv.getHeight2()-conv.getHeight1());
                double t2 = conv.getHeight1()
                        + Randomizer.nextDouble()*(conv.getHeight2()-conv.getHeight1());

                double tLower, tUpper;
                if (t1 < t2) {
                        tLower = t1;
                        tUpper = t2;
                } else {
                        tLower = t2;
                        tUpper = t1;
                }

                logHGF -= Math.log(0.5) - 2.0*Math.log(conv.getHeight2()-conv.getHeight1());

                // Select node below detour edge at random
                Node detour = acg.getNode(Randomizer.nextInt(acg.getNodeCount()-1));
                logHGF -= Math.log(1.0/(acg.getNodeCount()-1));

                // Abort if selected detour edge does not contain tLower and tUpper
                if (detour.getHeight() > tLower || detour.getParent().getHeight() < tUpper)
                        return Double.NEGATIVE_INFINITY;

                Conversion convA = new Conversion();
                convA.setNode1(conv.getNode1());
                convA.setHeight1(conv.getHeight1());
                convA.setNode2(detour);
                convA.setHeight2(tLower);
                logHGF -= drawAffectedRegion(convA);

                Conversion convB = new Conversion();
                convB.setNode1(detour);
                convB.setHeight1(tUpper);
                convB.setNode2(conv.getNode2());
                convB.setHeight2(conv.getHeight1());
                logHGF -= drawAffectedRegion(convB);

                acg.deleteConversion(conv);
                acg.addConversion(convA);
                acg.addConversion(convB);

                return logHGF;
        }

        /**
         * Detour deletion variant of move.
         * @return log HGF
         */
        private double removeDetour() {
                double logHGF = 0.0;

                return logHGF;
        }
}
