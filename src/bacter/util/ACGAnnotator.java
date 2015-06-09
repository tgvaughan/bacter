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

package bacter.util;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ACGAnnotator {

    public enum HeightStrategy { MEAN, MEDIAN }

    static class ACGAnnotatorOptions {
        public File inFile;
        public File outFile = new File("summary.tree");
        public double burninPercentage = 10.0;
        public HeightStrategy heightStrategy = HeightStrategy.MEAN;
    }

    public ACGAnnotator(ACGAnnotatorOptions options) {
    }

    public static ACGAnnotatorOptions getOptionsGUI() {

        ACGAnnotatorOptions options = new ACGAnnotatorOptions();
        boolean[] canceled = {false};

        JDialog dialog = new JDialog((JDialog)null, true);
        dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        dialog.setLocationRelativeTo(null);
        dialog.setTitle("ACGAnnotator");

        JLabel logFileLabel = new JLabel("ACG log file:");
        JLabel burninLabel = new JLabel("Burn-in percentage:");
        JLabel heightMethodLabel = new JLabel("Node height method:");

        JTextField inFilename = new JTextField(20);
        inFilename.setEditable(false);
        JButton infileButton = new JButton("Choose File");

        JSlider burninSlider = new JSlider(JSlider.HORIZONTAL,
                0, 100, 10);
        burninSlider.setMajorTickSpacing(50);
        burninSlider.setMinorTickSpacing(10);
        burninSlider.setPaintTicks(true);
        burninSlider.setPaintLabels(true);

        JComboBox<HeightStrategy> heightMethodCombo = new JComboBox<>(HeightStrategy.values());

        Container cp = dialog.getContentPane();

        JPanel mainPanel = new JPanel();

        GroupLayout layout = new GroupLayout(mainPanel);
        mainPanel.setLayout(layout);

        layout.setHorizontalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(burninLabel)
                        .addComponent(heightMethodLabel))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFilename)
                        .addComponent(burninSlider)
                        .addComponent(heightMethodCombo))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(infileButton)));

        layout.setVerticalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(inFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(infileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(burninLabel)
                        .addComponent(burninSlider,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE))
                .addGroup(layout.createParallelGroup()
                        .addComponent(heightMethodLabel)
                        .addComponent(heightMethodCombo,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)));

        cp.add(mainPanel, BorderLayout.CENTER);

        JPanel buttonPanel = new JPanel();

        JButton runButton = new JButton("Run");
        runButton.addActionListener((e) -> {
            options.burninPercentage = burninSlider.getValue();
            options.heightStrategy = (HeightStrategy)heightMethodCombo.getSelectedItem();
            dialog.setVisible(false);
        });
        runButton.setEnabled(false);
        buttonPanel.add(runButton);

        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener((e) -> {
            dialog.setVisible(false);
            canceled[0] = true;
        });
        buttonPanel.add(cancelButton);

        infileButton.addActionListener(e -> {
            JFileChooser chooser = new JFileChooser();
            chooser.setDialogTitle("Select ACG log file to summarize");
            int returnVal = chooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.inFile = chooser.getSelectedFile();
                inFilename.setText(chooser.getSelectedFile().getName());
                runButton.setEnabled(true);
            }
        });

        cp.add(buttonPanel, BorderLayout.PAGE_END);

                dialog.pack();
        dialog.setResizable(false);
        dialog.setVisible(true);

        if (canceled[0])
            return null;
        else
            return options;
    }

    public static void setupGUIOutput() {
        JFrame frame = new JFrame();
        frame.setTitle("ACGAnnotator");
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

        JTextArea textArea = new JTextArea(25, 80);
        textArea.setFont(new Font("monospaced", Font.PLAIN, 12));
        textArea.setEditable(false);
        frame.getContentPane().add(new JScrollPane(textArea), BorderLayout.CENTER);

        JButton closeButton = new JButton("Close");
        closeButton.addActionListener(e -> System.exit(0));
        JPanel buttonPanel = new JPanel();
        buttonPanel.add(closeButton);
        frame.getContentPane().add(buttonPanel, BorderLayout.PAGE_END);

        // Redirect streams to output window:
        OutputStream out = new OutputStream() {
            @Override
            public void write(int b) throws IOException {
                SwingUtilities.invokeLater(
                        () -> textArea.append(String.valueOf((char)b)));
            }
        };

        System.setOut(new PrintStream(out, true));
        System.setErr(new PrintStream(out, true));

        frame.pack();
        frame.setVisible(true);
    }

    public static String helpMessage =
            "ACGAnnotator - produces summaries of Bacter ACG log files.\n"
                    + "\n"
                    + "Usage: appstore ACGAnnotator [-help | [options] logFile [outputFile]\n"
                    + "\n"
                    + "Option                   Description\n"
                    + "--------------------------------------------------------------\n"
                    + "-help                    Display usage info.\n"
                    + "-heights {mean,median}   Choose node height method.\n"
                    + "-burnin percentage       Choose _percentage_ of log to discard\n"
                    + "                         in order to remove burn-in period.";

    public static void printUsageAndExit() {
        System.out.println(helpMessage);
        System.exit(0);
    }

    public static void printUsageAndError() {
        System.err.println("Error processing command line parameters.\n");
        System.err.println(helpMessage);
        System.exit(1);
    }

    public static ACGAnnotatorOptions getCLIOptions(String[] args) {
        ACGAnnotatorOptions options = new ACGAnnotatorOptions();

        int i=0;
        while (args[i].startsWith("-")) {
            switch(args[i]) {
                case "-help":
                    printUsageAndExit();
                    break;

                case "-burnin":
                    if (args.length<=i+1)
                        printUsageAndError();

                    try {
                        options.burninPercentage = Double.parseDouble(args[i+1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError();
                    }

                    if (options.burninPercentage<0 || options.burninPercentage>100)
                        printUsageAndError();

                    i += 1;
                    break;

                case "-heights":
                    if (args.length<=i+1)
                        printUsageAndError();

                    if (args[i+1].toLowerCase().equals("mean")) {
                        options.heightStrategy = HeightStrategy.MEAN;

                        i += 1;
                        break;
                    }

                    if (args[i+1].toLowerCase().equals("median")) {
                        options.heightStrategy = HeightStrategy.MEDIAN;

                        i += 1;
                        break;
                    }

                    printUsageAndError();

                default:
                    printUsageAndError();
            }

            i += 1;
        }

        if (i >= args.length)
            printUsageAndError();
        else
            options.inFile = new File(args[i]);

        if (i+1<args.length)
            options.outFile = new File(args[i+1]);

        return options;
    }

    public static void main(String[] args) {


        if (args.length == 0) {
            // Retrieve options from GUI:
            SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                    ACGAnnotatorOptions options = getOptionsGUI();

                    if (options == null)
                        System.exit(0);

                    setupGUIOutput();

                    // Run ACGAnnotator
                    new ACGAnnotator(options);
                }
            });

        } else {

            // Run ACGAnnotator
            new ACGAnnotator(getCLIOptions(args));
        }

    }
}
