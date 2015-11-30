package bacter.model;

import bacter.*;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;

import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Approximate ACG likelihood using pairwise distances.")
public class ACGLikelihoodApprox extends Distribution {

    public Input<ConversionGraph> acgInput = new Input<>(
            "acg",
            "Conversion graph.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> substRateInput = new Input<>(
            "substitutionRate",
            "Substitution rate.",
            Input.Validate.REQUIRED);

    public Input<Alignment> alignmentInput = new Input<>(
            "alignment",
            "Alignment corresponding to particular locus.",
            Input.Validate.REQUIRED);

    public Input<Locus> locusInput = new Input<>(
            "locus",
            "Locus alignment is associated with.",
            Input.Validate.REQUIRED);

    int nPairs;
    int[][] cumulativeHD;
    Alignment alignment;
    ConversionGraph acg;
    Locus locus;

    public ACGLikelihoodApprox() { }

    @Override
    public void initAndValidate() throws Exception {

        alignment = alignmentInput.get();
        acg = acgInput.get();
        locus = locusInput.get();

        // Pre-compute pairwise distance tables

        int nLeaves = acg.getLeafNodeCount();
        nPairs = nLeaves*(nLeaves+1)/2;
        cumulativeHD = new int[nPairs][alignment.getSiteCount()];

        for (int site=0; site<alignment.getSiteCount(); site++) {
            int pair = 0;
            for (int tIdx1=1; tIdx1<nLeaves; tIdx1++) {
                for (int tIdx2=0; tIdx2<tIdx1; tIdx2++) {
                    if (site==0)
                        cumulativeHD[pair][site] = 0;
                    else
                        cumulativeHD[pair][site] = cumulativeHD[pair][site-1];

                    int patternIdx = alignment.getPatternIndex(site);

                    if (alignment.getPattern(tIdx1, patternIdx)
                            != alignment.getPattern(tIdx2, patternIdx))
                        cumulativeHD[pair][site] += 1;

                    pair += 1;
                }
            }
        }
    }

    @Override
    public double calculateLogP() throws Exception {
        logP = 0.0;

        return logP;
    }

    double[] distancesTrue, distancesExpected;

    public void computeTrueDistances(Set<Conversion> activeConvs) {

    }

    public void computeExpectedDistances() {
        for (int pair=0; pair<nPairs; pair++) {

        }
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Sampling from this " +
                "ACGLikelihoodApprox is not supported.");
    }
}
