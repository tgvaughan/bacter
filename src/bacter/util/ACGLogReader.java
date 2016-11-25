package bacter.util;

import bacter.ConversionGraph;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public interface ACGLogReader extends Iterable<ConversionGraph> {

    int getACGCount();
    int getCorrectedACGCount();
}
