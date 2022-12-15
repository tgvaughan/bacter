package bacter.devutils;

import bacter.*;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;

import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Approximate ACG likelihood using pairwise distances.")
public class ComputeCPDMatrix extends beast.base.inference.Runnable {

    public Input<Alignment> alignmentInput = new Input<>(
            "alignment",
            "Alignment corresponding to particular locus.",
            Input.Validate.REQUIRED);

    public Input<String> outFileNameInput = new Input<>(
            "outFileName",
            "Name of output file",
            Input.Validate.REQUIRED);

    private int[][] cumulativeHD;
    private Alignment alignment;

    public ComputeCPDMatrix() { }

    @Override
    public void initAndValidate() {

        alignment = alignmentInput.get();
        computePairwiseDistances();
    }

    @Override
    public void run() throws Exception {

        try (PrintStream ps = new PrintStream(outFileNameInput.get())) {

            ps.print("site");

            for (int p=0; p<cumulativeHD.length; p++)
                ps.print("\tp" + p);

            ps.println();

            for (int site=0; site<alignment.getSiteCount(); site++) {

                ps.print(site);

                for (int p=0; p<cumulativeHD.length; p++) {

                        ps.print("\t" + cumulativeHD[p][site]);
                }
                ps.println();
            }

        } catch (IOException e) {
            System.err.println("Error writing output file '"
                    + outFileNameInput.get() + "'.");
        }
    }

    private void computePairwiseDistances() {
        // Pre-compute pairwise distance tables

        int nSeqs = alignment.getTaxonCount();
        int nPairs = nSeqs*(nSeqs-1)/2;
        cumulativeHD = new int[nPairs][alignment.getSiteCount()+1];

        for (int siteBoundary=0; siteBoundary<alignment.getSiteCount()+1; siteBoundary++) {
            int pair = 0;
            for (int tIdx1=0; tIdx1<nSeqs; tIdx1++) {
                for (int tIdx2=tIdx1+1; tIdx2<nSeqs; tIdx2++) {
                    if (siteBoundary == 0)
                        cumulativeHD[pair][siteBoundary] = 0;
                    else
                        cumulativeHD[pair][siteBoundary] = cumulativeHD[pair][siteBoundary - 1];

                    if (siteBoundary > 0) {
                        int patternIdx = alignment.getPatternIndex(siteBoundary - 1);

                        if (alignment.getPattern(tIdx1, patternIdx)
                                != alignment.getPattern(tIdx2, patternIdx))
                            cumulativeHD[pair][siteBoundary] += 1;
                    }

                    pair += 1;
                }
            }
        }
    }

}
