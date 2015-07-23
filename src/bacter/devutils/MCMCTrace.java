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

package bacter.devutils;

import beast.core.*;
import beast.core.util.CompoundDistribution;
import beast.core.util.Evaluator;
import beast.core.util.Log;
import beast.util.Randomizer;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("An extension of MCMC that generates a detailed trace of the " +
        "chain, including operators used, states proposed, Hastings ratios " +
        "and acceptance probabilities.")
public class MCMCTrace extends MCMC {

    public Input<String> stateTraceFileInput = new Input<>("stateTraceFile",
            "Name of file to which state trace will be written",
            Input.Validate.REQUIRED);

    public Input<List<Operator>> stateTraceOperatorsInput = new Input<>(
            "stateTraceOperator",
            "Operator to include in state trace.  If none are given, all" +
                    "operators are included.",
            new ArrayList<>());

    public Input<Integer> stateTraceStartInput = new Input<>("stateTraceStart",
            "Iteration at which state trace should begin.", 0);

    PrintStream stateTrace;
    List<Operator> stateTraceOperators;
    int stateTraceStart;

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        stateTrace = new PrintStream(stateTraceFileInput.get());
        stateTraceOperators = stateTraceOperatorsInput.get();
        stateTraceStart = stateTraceStartInput.get();
    }

    @Override
    protected void doLoop() throws Exception {
        int corrections = 0;
        if (burnIn > 0) {
            Log.warning.println("Please wait while BEAST takes " + burnIn + " pre-burnin samples");
        }

        List<String> targetDensityNames = new ArrayList<>();
        collectTargetDensityNames(posterior, targetDensityNames);

        List<Double> oldTargetDensities = new ArrayList<>();
        List<Double> newTargetDensities = new ArrayList<>();

        for (int sampleNr = -burnIn; sampleNr <= chainLength; sampleNr++) {
            final int currentState = sampleNr;


            state.store(currentState);

            final Operator operator = operatorSchedule.selectOperator();

            boolean stateTraceActive = (stateTraceOperators.isEmpty()
                    || stateTraceOperators.contains(operator))
                    && sampleNr >= stateTraceStart;

            if (stateTraceActive) {
                stateTrace.println("Iteration: " + sampleNr);

                stateTrace.print("Initial state:\n" + state.toXML(currentState));

                stateTrace.println("Operator: " + operator.getName());
            }

            final Distribution evaluatorDistribution = operator.getEvaluatorDistribution();
            Evaluator evaluator = null;

            if (evaluatorDistribution != null) {
                evaluator = () -> {
                    double logP = 0.0;

                    state.storeCalculationNodes();
                    state.checkCalculationNodesDirtiness();

                    try {
                        logP = evaluatorDistribution.calculateLogP();
                    } catch (Exception e) {
                        e.printStackTrace();
                        System.exit(1);
                    }

                    state.restore();
                    state.store(currentState);

                    return logP;
                };
            }

            if (stateTraceActive) {
                oldTargetDensities.clear();
                collectTargetDensities(posterior, oldTargetDensities);
            }

            final double logHastingsRatio = operator.proposal(evaluator);

            if (stateTraceActive)
                stateTrace.print("Proposed state: " + state.toXML(currentState));

            if (logHastingsRatio != Double.NEGATIVE_INFINITY) {

                if (operator.requiresStateInitialisation()) {
                    state.storeCalculationNodes();
                    state.checkCalculationNodesDirtiness();
                }

                newLogLikelihood = posterior.calculateLogP();

                logAlpha = newLogLikelihood - oldLogLikelihood + logHastingsRatio;

                if (stateTraceActive) {
//                    stateTrace.println("Target log density ratio: " + (newLogLikelihood - oldLogLikelihood));
                    newTargetDensities.clear();
                    collectTargetDensities(posterior, newTargetDensities);

                    stateTrace.print("Target log density ratios:");
                    for (int i=0; i<targetDensityNames.size(); i++) {
                        stateTrace.print(" " + targetDensityNames.get(i)
                        + ": " + (newTargetDensities.get(i)-oldTargetDensities.get(i)));
                    }
                    stateTrace.println();

                    stateTrace.println("Hastings Ratio: " + logHastingsRatio);
                    stateTrace.println("Log acceptance probability: " + Math.min(logAlpha, 0.0));
                }
                if (logAlpha >= 0 || Randomizer.nextDouble() < Math.exp(logAlpha)) {
                    // accept
                    oldLogLikelihood = newLogLikelihood;
                    state.acceptCalculationNodes();

                    if (sampleNr >= 0) {
                        operator.accept();
                    }

                    if (stateTraceActive)
                        stateTrace.println("Result: ACCEPT");
                } else {
                    // reject
                    if (sampleNr >= 0) {
                        operator.reject(newLogLikelihood == Double.NEGATIVE_INFINITY ? -1 : 0);
                    }
                    state.restore();
                    state.restoreCalculationNodes();

                    if (stateTraceActive)
                        stateTrace.println("Result: REJECT");
                }
                state.setEverythingDirty(false);
            } else {
                // operation failed
                if (sampleNr >= 0) {
                    operator.reject(-2);
                }
                state.restore();
                if (!operator.requiresStateInitialisation()) {
                    state.setEverythingDirty(false);
                    state.restoreCalculationNodes();
                }

                if (stateTraceActive)
                    stateTrace.println("Result: DIRECT REJECT");
            }
            log(sampleNr);

            if (sampleNr >= 0)
                operator.optimize(logAlpha);

            // make sure we always save just before exiting
            if (storeEvery > 0 && (sampleNr + 1) % storeEvery == 0 || sampleNr == chainLength) {
                /*final double fLogLikelihood = */
                state.robustlyCalcNonStochasticPosterior(posterior);
                state.storeToFile(sampleNr);
                operatorSchedule.storeToFile();
            }

            if (stateTraceActive)
                stateTrace.println("\n===\n");
        }
        if (corrections > 0) {
            System.err.println("\n\nNB: " + corrections + " posterior calculation corrections were required. This analysis may not be valid!\n\n");
        }
    }

    protected void collectTargetDensities(final Distribution distr, List<Double> densities) {
        densities.add(distr.getCurrentLogP());
        if (distr instanceof CompoundDistribution) {
            for (final Distribution childDistr : ((CompoundDistribution) distr).pDistributions.get()) {
                collectTargetDensities(childDistr, densities);
            }
        }
    }

    protected void collectTargetDensityNames(final Distribution distr, List<String> densityNames) {
        densityNames.add(distr.getID());
        if (distr instanceof CompoundDistribution) {
            for (final Distribution childDistr : ((CompoundDistribution) distr).pDistributions.get()) {
                collectTargetDensityNames(childDistr, densityNames);
            }
        }
    }
}
