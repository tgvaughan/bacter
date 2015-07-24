package bacter.devutils;

import beast.core.*;
import beast.core.Runnable;
import beast.util.XMLParser;

import java.io.File;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class OperatorTester {

    public static void main(String[] args) {

        if (args.length < 4) {
            System.out.println("Usage: OperatorTester model.xml operator_ID state_file iterations");
            System.exit(0);
        }

        String modelFileName = args[0];
        String operatorID = args[1];
        String stateFileName = args[2];
        int Niter = Integer.parseInt(args[3]);

        XMLParser xmlParser = new XMLParser();
        Runnable runnable = null;
        try {
            runnable = xmlParser.parseFile(new File(modelFileName));
        } catch (Exception e) {
            System.out.println("Error parsing model XML file.");
            e.printStackTrace();
        }

        if (runnable == null)
            System.exit(1);

        if (!(runnable instanceof MCMC)) {
            System.out.println("XML file does not describe MCMC analysis.");
            System.exit(1);
        }

        MCMC mcmc = (MCMC)runnable;

        Operator operator = null;
        for (Operator thisOperator : mcmc.operatorsInput.get()) {
            if (thisOperator.getID().equals(operatorID)) {
                operator = thisOperator;
                break;
            }
        }

        if (operator == null) {
            System.out.println("Model does not include operator with ID "
                    + operatorID + ".");
            System.exit(1);
        }

        State state = mcmc.startStateInput.get();

        state.setStateFileName(stateFileName);

        try {
            state.restoreFromFile();
        } catch (Exception e) {
            System.out.println("Error loading state from file.");
            e.printStackTrace();
            System.exit(1);
        }

        try {
            state.robustlyCalcPosterior(mcmc.posteriorInput.get());
        } catch (Exception e) {
            System.out.println("Failed to calculate posterior.");
            System.exit(1);
        }

        Logger.FILE_MODE = Logger.LogFileMode.overwrite;

        for (Logger logger : mcmc.loggersInput.get()) {
            try {
                logger.everyInput.setValue(1, logger);
                logger.initAndValidate();
                logger.init();
                logger.log(0);
            } catch (Exception e) {
                System.out.println("Error initializing logger.");
                System.exit(1);
            }
        }

        for (int i=1; i<Niter; i++) {
            try {
                state.robustlyCalcPosterior(mcmc.posteriorInput.get());
            } catch (Exception e) {
                System.out.println("Failed to calculate posterior.");
                System.exit(1);
            }

            operator.proposal();

            state.storeCalculationNodes();
            state.checkCalculationNodesDirtiness();

            try {
                mcmc.posteriorInput.get().calculateLogP();
            } catch (Exception e) {
                System.out.println("Failed to calculate posterior.");
                System.exit(1);
            }

            for (Logger logger : mcmc.loggersInput.get())
                logger.log(i);

            operator.reject();
            state.restore();
            state.restoreCalculationNodes();
        }

        for (Logger logger : mcmc.loggersInput.get())
            logger.close();
    }
}
