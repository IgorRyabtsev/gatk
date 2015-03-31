package org.broadinstitute.hellbender.utils.sam.markduplicates;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SamRecordWithOrdinal;

/**
 * This class sets the duplicate read flag as the result state when examining sets of records.
 *
 * @author nhomer
 */
public class SamRecordWithOrdinalAndSetDuplicateReadFlag extends SamRecordWithOrdinal {

    public SamRecordWithOrdinalAndSetDuplicateReadFlag() {
        super();
    }

    public SamRecordWithOrdinalAndSetDuplicateReadFlag(final SAMRecord record, final long recordIndex) {
        super(record, recordIndex);
    }

    @Override
    public void setResultState(final boolean resultState) {
        this.getRecord().setDuplicateReadFlag(resultState);
    }
}