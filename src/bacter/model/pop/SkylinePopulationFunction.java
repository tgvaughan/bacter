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

import java.util.*;

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

    private boolean dirty;

    double[] intensities;
    double[] groupBoundaries;

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

        dirty = true;
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
    public void prepare() {

        if (!dirty)
            return;

        List<Event> cfEvents = acg.getCFEvents();
        groupBoundaries = new double[groupSizes.getDimension() - 1];
        int lastEventIdx = -1;
        for (int i=0; i<groupBoundaries.length; i++) {
            int cumulant = 0;
            do {
                lastEventIdx += 1;
                if (cfEvents.get(lastEventIdx).getType() == CFEventList.EventType.COALESCENCE)
                    cumulant += 1;
            } while (cumulant < groupSizes.getValue(i));

            groupBoundaries[i] = cfEvents.get(lastEventIdx).getHeight();
        }

        intensities = new double[groupSizes.getDimension()+1];
        intensities[0] = 0.0;
        double lastBoundary = 0.0;

        for (int i=1; i<intensities.length; i++) {
            intensities[i] = intensities[i-1]
                    + (groupBoundaries[i-1]-lastBoundary)/popSizes.getValue(i-1);
        }

        dirty = false;
    }

    @Override
    protected boolean requiresRecalculation() {
        dirty = true;
        return true;
    }

    @Override
    protected void store() {
        dirty = true;
        super.store();
    }

    @Override
    protected void restore() {
        dirty = true;
        super.restore();
    }

    @Override
    public double getPopSize(double t) {
        prepare();

        if (t < 0)
            return getPopSize(0);

        int interval = Arrays.binarySearch(groupBoundaries, t);

        if (interval<0)
            interval = -(interval + 1) - 1;  // boundary to the left of time.

        return popSizes.getValue(interval+1);
    }

    @Override
    public double getIntensity(double t) {
        prepare();

        if (t < 0 )
            return -t/popSizes.getValue(0);

        int interval = Arrays.binarySearch(groupBoundaries, t);

        if (interval<0)
            interval = -(interval + 1) - 1; // boundary to the left of time.

        double intensityGridTime = interval >= 0
                ? groupBoundaries[interval]
                : 0.0;

        return intensities[interval+1] + (t-intensityGridTime)/popSizes.getValue(interval+1);
    }

    @Override
    public double getInverseIntensity(double x) {
        throw new RuntimeException("Method not implemented");
    }
}
