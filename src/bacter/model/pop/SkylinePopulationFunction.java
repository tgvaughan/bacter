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

package bacter.model.pop;

import bacter.CFEventList;
import bacter.CFEventList.Event;
import bacter.ConversionGraph;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.PopulationFunction;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * In BEAST 2, BSP is implemented as a tree distribution rather than
 * a population function.  The population function approach is much
 * more flexible and is directly applicable to BACTER's model.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SkylinePopulationFunction extends PopulationFunction.Abstract {

    public Input<ConversionGraph> acgInput = new Input<>("acg",
            "Conversion graph", Input.Validate.REQUIRED);

    public Input<RealParameter> popSizesInput = new Input<>("popSizes",
            "Population sizes in intervals", Input.Validate.REQUIRED);

    public Input<IntegerParameter> groupSizesInput = new Input<>("groupSizes",
            "Group size parameter.", Input.Validate.REQUIRED);

    public Input<Boolean> initGroupSizesInput = new Input<>("initGroupSizes",
            "If true (default), initialize group sizes parameter.", true);

    ConversionGraph acg;
    RealParameter popSizes;
    IntegerParameter groupSizes;

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        acg = acgInput.get();
        popSizes = popSizesInput.get();
        groupSizes = groupSizesInput.get();

        // Initialize groupSizes to something sensible.
        if (initGroupSizesInput.get()) {
            int nCFEvents = acg.getCFEvents().size();
            Integer[] values = new Integer[groupSizes.getDimension()];
            int cumulant = 0;
            for (int i=0; i<values.length; i++) {
                values[i] = nCFEvents/values.length;
                cumulant += values[i];
            }
            values[values.length-1] += nCFEvents - cumulant;

            IntegerParameter newParam = new IntegerParameter(values);
            groupSizes.assignFromWithoutID(newParam);
        }
    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        ids.add(acgInput.get().getID());
        ids.add(popSizesInput.get().getID());
        ids.add(groupSizesInput.get().getID());

        return ids;
    }

    @Override
    public double getPopSize(double t) {
        List<Event> cfEvents = acg.getCFEvents();

        int interval = Collections.binarySearch(cfEvents, new Event(t),
                new Comparator<Event>() {
                    @Override
                    public int compare(Event e1, Event e2) {
                        if (e1.getHeight()<e2.getHeight())
                            return -1;

                        if (e1.getHeight()>e2.getHeight())
                            return 1;

                        return 0;
                    }
                });

        if (interval<0)
            interval = -(interval + 1) - 1;  // event to the left of time.

        if (interval == -1)
            return popSizes.getValue(0);

        int cumulant = 0;
        for (int i=0; i<groupSizes.getDimension(); i++) {
            cumulant += groupSizes.getValue(i);
            if (cumulant>interval)
                return popSizes.getValue(i);
        }

        return popSizes.getValue(popSizes.getDimension()-1);
    }

    @Override
    public double getIntensity(double t) {
        return 0;
    }

    @Override
    public double getInverseIntensity(double x) {
        throw new RuntimeException("Method not implemented");
    }
}
