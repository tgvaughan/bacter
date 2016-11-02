/*
 * Copyright (C) 2015 Tim Vaughan (tgvaughan@gmail.com)
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
package bacter.operators;

import bacter.Conversion;
import bacter.Locus;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;

/**
 * Abstract class of ACG operators that use the clonal origin model as the
 * basis for adding new converted edges and their affected sites to an
 * existing ConversionGraph.
 *
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
public abstract class ConversionCreationOperator extends EdgeCreationOperator {

    public Input<RealParameter> deltaInput = new Input<>(
            "delta",
            "Tract length parameter.",
            Input.Validate.REQUIRED);

    double[] relativeLocusSizes;

    @Override
    public void initAndValidate() {

        relativeLocusSizes = new double[acgInput.get().getLoci().size()];
        for (int i=0; i<acgInput.get().getLoci().size(); i++)
            relativeLocusSizes[i] = acgInput.get().getLoci().get(i).getSiteCount()/
                    (double)acgInput.get().getTotalSequenceLength();

        super.initAndValidate();
    }


    /**
     * Choose region to be affected by this conversion.
     *
     * @param conv Conversion object whose region is to be set.
     * @return log probability density of chosen attachment.
     */
    public double drawAffectedRegion(Conversion conv) {
        if (!acg.wholeLocusModeOn())
            return drawAffectedRegionUnrestricted(conv);
        else
            return drawAffectedRegionRestricted(conv);
    }

    /**
     * Calculate probability of choosing region affected by the
     * given conversion.
     *
     * @param conv conversion region is associated with
     * @return log probability density
     */
    public double getAffectedRegionProb(Conversion conv) {
        if (!acg.wholeLocusModeOn())
            return getAffectedRegionProbUnrestricted(conv);
        else
            return getAffectedRegionProbRestricted(conv);
    }

    /**
     * Choose region to be affected by this conversion using the
     * full ClonalOrigin model.
     * 
     * @param conv Conversion object where these sites are stored.
     * @return log probability density of chosen attachment.
     */
    protected double drawAffectedRegionUnrestricted(Conversion conv) {
        double logP = 0.0;

        // Total effective number of possible start sites
        double alpha = acg.getTotalSequenceLength()
                + acg.getLoci().size()*(deltaInput.get().getValue() - 1.0);

        // Draw location of converted region.
        int startSite = -1;
        int endSite;
        Locus locus = null;

        double u = Randomizer.nextDouble()*alpha;
        for (Locus thisLocus : acg.getLoci()) {
            if (u < deltaInput.get().getValue() - 1.0 + thisLocus.getSiteCount()) {
                locus = thisLocus;

                if (u < deltaInput.get().getValue()) {
                    startSite = 0;
                    logP += Math.log(deltaInput.get().getValue() / alpha);
                } else {
                    startSite = (int)Math.ceil(u - deltaInput.get().getValue());
                    logP += Math.log(1.0 / alpha);
                }

                break;
            }

            u -= deltaInput.get().getValue() - 1.0 + thisLocus.getSiteCount();
        }

        if (locus == null)
            throw new IllegalStateException("Programmer error: " +
                    "loop in drawAffectedRegion() fell through.");

        endSite = startSite + (int)Randomizer.nextGeometric(1.0/deltaInput.get().getValue());
        endSite = Math.min(endSite, locus.getSiteCount()-1);

        // Probability of end site:
        if (endSite == locus.getSiteCount()-1) {
            logP += (locus.getSiteCount()-1-startSite)
                    *Math.log(1.0 - 1.0/deltaInput.get().getValue());
        } else {
            logP += (endSite-startSite)
                    *Math.log(1.0 - 1.0/deltaInput.get().getValue())
                    -Math.log(deltaInput.get().getValue());
        }

        conv.setLocus(locus);
        conv.setStartSite(startSite);
        conv.setEndSite(endSite);

        return logP;
    }
    
    /**
     * Calculate probability of choosing region affected by the given
     * conversion under the full ClonalOrigin model.
     * 
     * @param conv conversion region is associated with
     * @return log probability density
     */
    protected double getAffectedRegionProbUnrestricted(Conversion conv) {
        double logP = 0.0;

        // Total effective number of possible start sites
        double alpha = acg.getTotalSequenceLength()
                + acg.getLoci().size()*(deltaInput.get().getValue() - 1.0);

        // Calculate probability of converted region.
        if (conv.getStartSite()==0)
            logP += Math.log(deltaInput.get().getValue() / alpha);
        else
            logP += Math.log(1.0 / alpha);

        // Probability of end site:
        if (conv.getEndSite() == conv.getLocus().getSiteCount()-1) {
            logP += (conv.getLocus().getSiteCount()-1-conv.getStartSite())
                    *Math.log(1.0 - 1.0/deltaInput.get().getValue());
        } else {
            logP += (conv.getEndSite()-conv.getStartSite())
                    *Math.log(1.0 - 1.0/deltaInput.get().getValue())
                    -Math.log(deltaInput.get().getValue());
        }

        return logP;
    }

    /**
     * Select affected region when region endpoint restriction is in place.
     *
     * @param conv Conversion object where these sites are stored.
     * @return log probability density of chosen attachment.
     */
    protected double drawAffectedRegionRestricted(Conversion conv) {
        int locusIdx = Randomizer.randomChoicePDF(relativeLocusSizes);
        Locus locus = acg.getLoci().get(locusIdx);

        conv.setLocus(locus);
        conv.setStartSite(0);
        conv.setEndSite(locus.getSiteCount()-1);

        return Math.log(relativeLocusSizes[locusIdx]);
    }

    /**
     * Calculate probability of choosing region affected by the given
     * conversion when region endpoint restriction is in place.
     *
     * @param conv conversion region is associated with
     * @return log probability density
     */
    protected double getAffectedRegionProbRestricted(Conversion conv) {
        if (conv.getStartSite()>0 || conv.getEndSite()<conv.getLocus().getSiteCount()-1)
            return Double.NEGATIVE_INFINITY;

        return Math.log(relativeLocusSizes[acg.getLoci().indexOf(conv.getLocus())]);
    }
}
