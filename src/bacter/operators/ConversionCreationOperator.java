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
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import beast.base.util.GammaFunction;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.random.MersenneTwister;

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

        relativeLocusSizes = new double[acgInput.get().getConvertibleLoci().size()];
        for (int i = 0; i<acgInput.get().getConvertibleLoci().size(); i++)
            relativeLocusSizes[i] = acgInput.get().getConvertibleLoci().get(i).getSiteCount()/
                    (double)acgInput.get().getTotalConvertibleSequenceLength();

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
        //circular genome mode edit
        double alpha = acg.getTotalConvertibleSequenceLength()                      
                + (acg.circularGenomeModeOn() ? 0 : acg.getConvertibleLoci().size()*(deltaInput.get().getValue() - 1.0));

        // Draw location of converted region.
        //circular genome mode edit
        int startSite = -1;
        int endSite;
        Locus locus = null;

        double u = Randomizer.nextDouble()*alpha;

        if (!acg.circularGenomeModeOn()) {
            for (Locus thisLocus : acg.getConvertibleLoci()) {
                if (u < deltaInput.get().getValue() - 1.0 + thisLocus.getSiteCount()) {
                    locus = thisLocus;
                    if (u < deltaInput.get().getValue()) {
                        startSite = 0;
                        logP += Math.log(deltaInput.get().getValue() / alpha);
                    } else {
                        startSite = (int) Math.ceil(u - deltaInput.get().getValue());
                        logP += Math.log(1.0 / alpha);
                    }

                    break;
                }

                u -= deltaInput.get().getValue() - 1.0 + thisLocus.getSiteCount();
            }
        } else {
            locus = acg.getConvertibleLoci().get(0); //only one locus if genome circular
            startSite = (int) u;
            logP += Math.log(1.0 / alpha);
        }

        if (locus == null)
            throw new IllegalStateException("Programmer error: " +
                    "loop in drawAffectedRegion() fell through.");

		//circular genome mode edit
        if (!acg.circularGenomeModeOn()) {
            endSite = startSite + (int) Randomizer.nextGeometric(1.0 / deltaInput.get().getValue());
            endSite = Math.min(endSite, locus.getSiteCount() - 1);
        } else if (!acg.endSiteBetaBinomOn()) {
            int convLength;
            convLength = acg.getTotalConvertibleSequenceLength();
            while (convLength >= 0.5*acg.getTotalConvertibleSequenceLength()){
                convLength = (int) Randomizer.nextGeometric(1.0 / deltaInput.get().getValue());
            }
            endSite = ((startSite + convLength) >= acg.getTotalConvertibleSequenceLength()) ?
                    (startSite - acg.getTotalConvertibleSequenceLength() + convLength) : (startSite + convLength);

        } else {
            MersenneTwister rng = new MersenneTwister();
            int numTrials = (int) Math.floor((acg.getTotalConvertibleSequenceLength() - 1.) * 0.5);
            rng.setSeed(Randomizer.nextInt());
            BetaDistribution beta_dist = new BetaDistribution(rng, numTrials/(numTrials-deltaInput.get().getValue()),
                    numTrials/deltaInput.get().getValue(), 1.0E-9D);
            int convLength;
            BinomialDistribution binom_dist = new BinomialDistribution(rng, numTrials, beta_dist.sample());
            convLength = binom_dist.sample();
            endSite = ((startSite + convLength) >= acg.getTotalConvertibleSequenceLength()) ?
                    (startSite - acg.getTotalConvertibleSequenceLength() + convLength) : (startSite + convLength);
        }

        // Probability of end site:
        //circular genome mode edit
        if (acg.circularGenomeModeOn()) {
            if (acg.endSiteBetaBinomOn()) {
                int halfGenomeLength = (int) Math.floor((acg.getTotalConvertibleSequenceLength() - 1.) * 0.5);
                int kBetaBinom = (startSite <= endSite) ? (endSite - startSite) : (acg.getTotalConvertibleSequenceLength() - startSite + endSite);
                double aBetaBinom = halfGenomeLength / (halfGenomeLength - deltaInput.get().getValue());
                double bBetaBinom = halfGenomeLength / deltaInput.get().getValue();
                logP += GammaFunction.lnGamma(halfGenomeLength + 1) - GammaFunction.lnGamma(kBetaBinom + 1) - GammaFunction.lnGamma(halfGenomeLength - kBetaBinom + 1)
                        + GammaFunction.lnGamma(kBetaBinom + aBetaBinom) + GammaFunction.lnGamma(halfGenomeLength - kBetaBinom + bBetaBinom) - GammaFunction.lnGamma(bBetaBinom)
                        - GammaFunction.lnGamma(halfGenomeLength + aBetaBinom + bBetaBinom) + GammaFunction.lnGamma(aBetaBinom + bBetaBinom) - GammaFunction.lnGamma(aBetaBinom);
            } else {
                logP += ((startSite <= endSite) ? (endSite - startSite) : (acg.getTotalConvertibleSequenceLength() - startSite + endSite))
                        *Math.log(1.0 - 1.0/deltaInput.get().getValue())
                        -Math.log(deltaInput.get().getValue())
                        -Math.log(1.0-Math.pow(1.0-1.0/deltaInput.get().getValue(), (int) Math.floor((acg.getTotalConvertibleSequenceLength()) * 0.5)));
            }
        } else if (endSite == locus.getSiteCount()-1) {
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
        //circular genome mode edit
        double alpha;
        alpha = acg.getTotalConvertibleSequenceLength()
                + ( acg.circularGenomeModeOn() ? 0 : acg.getConvertibleLoci().size() * (deltaInput.get().getValue() - 1.0) );

        // Calculate probability of converted region.
        //circular genome mode edit
        if (conv.getStartSite()==0 && !acg.circularGenomeModeOn())
            logP += Math.log(deltaInput.get().getValue() / alpha);
        else
            logP += Math.log(1.0 / alpha);

        // Probability of end site:
        //circular genome mode edit
        if (acg.circularGenomeModeOn()) {
            if (acg.endSiteBetaBinomOn()) {
                int halfGenomeLength = (int) Math.floor((acg.getTotalConvertibleSequenceLength() - 1.) * 0.5);
                int kBetaBinom = conv.getSiteCount() - 1;
                double aBetaBinom = halfGenomeLength / (halfGenomeLength - deltaInput.get().getValue());
                double bBetaBinom = halfGenomeLength / deltaInput.get().getValue();
                logP += GammaFunction.lnGamma(halfGenomeLength + 1) - GammaFunction.lnGamma(kBetaBinom + 1) - GammaFunction.lnGamma(halfGenomeLength - kBetaBinom + 1)
                        + GammaFunction.lnGamma(kBetaBinom + aBetaBinom) + GammaFunction.lnGamma(halfGenomeLength - kBetaBinom + bBetaBinom) - GammaFunction.lnGamma(bBetaBinom)
                        - GammaFunction.lnGamma(halfGenomeLength + aBetaBinom + bBetaBinom) + GammaFunction.lnGamma(aBetaBinom + bBetaBinom) - GammaFunction.lnGamma(aBetaBinom);
            } else {
                logP += (conv.getSiteCount() - 1)
                        *Math.log(1.0 - 1.0/deltaInput.get().getValue())
                        -Math.log(deltaInput.get().getValue())
                        -Math.log(1.0-Math.pow(1.0-1.0/deltaInput.get().getValue(), (int) Math.floor((acg.getTotalConvertibleSequenceLength()) * 0.5)));
            }
        } else if (conv.getEndSite() == conv.getLocus().getSiteCount()-1) {
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
        Locus locus = acg.getConvertibleLoci().get(locusIdx);

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

        return Math.log(relativeLocusSizes[acg.getConvertibleLoci().indexOf(conv.getLocus())]);
    }
}
