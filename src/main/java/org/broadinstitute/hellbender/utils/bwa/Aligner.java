package org.broadinstitute.hellbender.utils.bwa;

import java.util.List;

public interface Aligner extends AutoCloseable {
    // return the reference contig names
    List<String> getReferenceContigNames();

    // align each sequence, and return a list of Alignments for each
    List<List<Alignment>> alignSeqs( final List<byte[]> sequences );

    // an aligner is an AutoCloseable that doesn't throw checked exceptions
    @Override void close();
}
