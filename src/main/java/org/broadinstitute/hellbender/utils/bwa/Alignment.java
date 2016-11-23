package org.broadinstitute.hellbender.utils.bwa;

import htsjdk.samtools.*;
import htsjdk.samtools.util.QualityUtil;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.Arrays;
import java.util.List;

public class Alignment {
    private final int samFlag;     // flag bits per SAM format standard
    private final int refId;       // index into reference dictionary (-1 if unmapped)
    private final int refStart;    // 0-based coordinate, inclusive (-1 if unmapped)
    private final int refEnd;      // 0-based coordinate, exclusive (-1 if unmapped)
    private final int tigStart;    // 0-based coordinate, inclusive (-1 if unmapped)
    private final int tigEnd;      // 0-based coordinate, exclusive (-1 if unmapped)
    private final int mapQual;     // phred-scaled mapping quality (-1 if unmapped)
    private final int nMismatches; // number of mismatches (i.e., value of the NM tag in a SAM/BAM) (-1 if unmapped)
    private final int alignerScore; // for AS tag
    private final int suboptimalScore; // for bwa-specific XS tag
    private final Cigar cigar;     // cigar for alignment (empty if unmapped)
    private final String mdTag;    // the MD tag
    private final String xaTag;    // the XA tag
    private final int mateRefId;   // mate's refId (-1 if unpaired or if mate unmapped)
    private final int mateRefStart;// mate's reference start (-1 if unpaired or if mate unmapped)
    private final int templateLen; // inferred template length (0 if unpaired, mate unmapped, or on different ref contigs)

    public Alignment( final int samFlag, final int refId, final int refStart, final int refEnd,
                      final int tigStart, final int tigEnd, final int mapQual,
                      final int nMismatches, final int alignerScore, final int suboptimalScore,
                      final Cigar cigar, final String mdTag, final String xaTag,
                      final int mateRefId, final int mateRefStart, final int templateLen ) {
        this.samFlag = samFlag;
        this.refId = refId;
        this.refStart = refStart;
        this.refEnd = refEnd;
        this.tigStart = tigStart;
        this.tigEnd = tigEnd;
        this.mapQual = mapQual;
        this.nMismatches = nMismatches;
        this.alignerScore = alignerScore;
        this.suboptimalScore = suboptimalScore;
        this.cigar = cigar;
        this.mdTag = mdTag;
        this.xaTag = xaTag;
        this.mateRefId = mateRefId;
        this.mateRefStart = mateRefStart;
        this.templateLen = templateLen;
    }

    public int getSamFlag() { return samFlag; }
    public int getRefId() { return refId; }
    public int getRefStart() { return refStart; }
    public int getRefEnd() { return refEnd; }
    public int getTigStart() { return tigStart; }
    public int getTigEnd() { return tigEnd; }
    public int getMapQual() { return mapQual; }
    public int getNMismatches() { return nMismatches; }
    public int getAlignerScore() { return alignerScore; }
    public int getSuboptimalScore() { return suboptimalScore; }
    public Cigar getCigar() { return cigar; }
    public String getMDTag() { return mdTag; }
    public int getMateRefId() { return mateRefId; }
    public int getMateRefStart() { return mateRefStart; }
    public int getTemplateLen() { return templateLen; }

    public GATKRead apply( final GATKRead originalRead, final List<String> refNames, final SAMFileHeader header,
                           final boolean softclipAlts, final boolean justLikeBWA ) {
        final SAMRecord samRecord;
        try {
            samRecord = (SAMRecord)originalRead.convertToSAMRecord(header).clone();
        } catch ( CloneNotSupportedException x ) {
            throw new GATKException("can't clone SAMRecord", x);
        }
        samRecord.setFlags(samFlag);
        if ( refId >= 0 ) samRecord.setReferenceName(refNames.get(refId));
        else if ( mateRefId >= 0 ) samRecord.setReferenceName(refNames.get(mateRefId));
        if ( refStart >= 0 ) samRecord.setAlignmentStart(refStart+1);
        else if ( mateRefStart >= 0 ) samRecord.setAlignmentStart(mateRefStart+1);
        if ( mapQual >= 0 ) samRecord.setMappingQuality(mapQual);
        byte[] seq = samRecord.getReadBases();
        byte[] quals = samRecord.getBaseQualities();
        if ( (samFlag & SAMFlag.READ_REVERSE_STRAND.intValue()) != 0 &&
                (samFlag & SAMFlag.NOT_PRIMARY_ALIGNMENT.intValue()) == 0 ) {
            seq = BaseUtils.simpleReverseComplement(seq);
            quals = Arrays.copyOf(quals, quals.length);
            SequenceUtil.reverseQualities(quals);
        }
        if ( cigar != null && !cigar.isEmpty() ) {
            Cigar tmpCigar = cigar;
            if ( !softclipAlts && (samFlag & SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue()) != 0 ) {
                if ( tmpCigar.getFirstCigarElement().getOperator() == CigarOperator.S ||
                        tmpCigar.getLastCigarElement().getOperator() == CigarOperator.S ) {
                    tmpCigar = new Cigar();
                    for ( final CigarElement ele : cigar ) {
                        if ( ele.getOperator() == CigarOperator.S ) {
                            tmpCigar.add(new CigarElement(ele.getLength(), CigarOperator.H));
                        }
                        else {
                            tmpCigar.add(ele);
                        }
                    }
                }
                seq = Arrays.copyOfRange(seq, tigStart, tigEnd);
                quals = Arrays.copyOfRange(quals, tigStart, tigEnd);
            }
            samRecord.setCigar(tmpCigar);
            samRecord.setAttribute("NM", nMismatches);
            samRecord.setAttribute("AS", alignerScore);
            samRecord.setAttribute("XS", suboptimalScore);
            samRecord.setAttribute("MD", mdTag);
            samRecord.setAttribute("XA", xaTag);
        }
        else if ( justLikeBWA ) {
            samRecord.setAttribute("AS", 0);
            samRecord.setAttribute("XS", 0);
        }
        if ( mateRefId >= 0 ) samRecord.setMateReferenceName(refNames.get(mateRefId));
        else if ( refId >= 0 ) samRecord.setMateReferenceName(refNames.get(refId));
        if ( mateRefStart >= 0 ) samRecord.setMateAlignmentStart(mateRefStart+1);
        else if ( refStart >= 0 ) samRecord.setMateAlignmentStart(refStart+1);
        if ( templateLen != 0 ) samRecord.setInferredInsertSize(templateLen);
        if ( (samFlag & SAMFlag.NOT_PRIMARY_ALIGNMENT.intValue()) == 0 ) {
            samRecord.setReadBases(seq);
            samRecord.setBaseQualities(quals);
        } else {
            samRecord.setReadBases(SAMRecord.NULL_SEQUENCE);
            samRecord.setBaseQualities(SAMRecord.NULL_QUALS);
        }
        return SAMRecordToGATKReadAdapter.headerlessReadAdapter(samRecord);
    }

    public String asTag( final List<String> refNames ) {
        return refNames.get(refId)+","+(refStart+1)+","+((samFlag&SAMFlag.READ_REVERSE_STRAND.intValue())!=0?"-":"+")+","+
                cigar.toString()+","+mapQual+","+nMismatches+";";
    }
}
