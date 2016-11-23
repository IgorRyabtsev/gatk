package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFlag;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.bwa.Aligner;
import org.broadinstitute.hellbender.utils.bwa.Alignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import scala.Tuple2;

import java.io.IOException;
import java.util.*;

import static org.broadinstitute.hellbender.tools.spark.sv.ContigsCollection.ContigID;
import static org.broadinstitute.hellbender.tools.spark.sv.ContigsCollection.ContigSequence;

public class ContigAligner {

    private final String indexImageFile;

    public ContigAligner(final String indexImageFile) {
        this.indexImageFile = indexImageFile;
    }

    /**
     * Takes a collection of assembled contigs and aligns them to the reference with jBWA. Non-canonical
     * (secondary) alignments are filtered out, preserving the primary and supplementary alignments.
     * Within the output list, alignments are sorted first by contig (based on the order in which
     * the contigs were passed in, and then by their start position on the contig).
     *
     * @param assemblyId An identifier for the assembly or set of contigs
     * @param contigsCollection The set of all canonical (primary or supplementary) alignments for the contigs.
     */
    public List<AlignmentRegion> alignContigs(final String assemblyId, final ContigsCollection contigsCollection) {
        final List<AlignmentRegion> alignedContigs = new ArrayList<>(contigsCollection.getContents().size());
        try ( final BwaMemIndex.BwaMemAligner aligner = BwaMemIndex.getInstance(indexImageFile).createAligner() ) {
            final List<String> refNames = aligner.getReferenceContigNames();
            final List<Tuple2<ContigID, ContigSequence>> contents = contigsCollection.getContents();
            final List<byte[]> seqs = new ArrayList<>(contents.size());
            for ( final Tuple2<ContigID, ContigSequence> contigInfo : contents ) {
                seqs.add(contigInfo._2().toString().getBytes());
            }
            Iterator<List<Alignment>> alignmentsItr = aligner.alignSeqs(seqs).iterator();
            for ( final Tuple2<ContigID, ContigSequence> contigInfo : contents ) {
                final String contigId = contigInfo._1.toString();
                final int contigLen = contigInfo._2().toString().length();
                final List<Alignment> alignments = alignmentsItr.next();

                // filter out secondary alignments, convert to AlignmentRegion objects and sort by alignment start pos
                alignments.stream()
                        .filter(a -> (a.getSamFlag()&SAMFlag.NOT_PRIMARY_ALIGNMENT.intValue())==0)
                        .filter(a -> (a.getSamFlag()&SAMFlag.READ_UNMAPPED.intValue())==0)
                        .map(a -> new AlignmentRegion(assemblyId, contigId, contigLen, a, refNames))
                        .sorted(Comparator.comparing(a -> a.startInAssembledContig))
                        .forEach(alignedContigs::add);
            }
        }

        return alignedContigs;
    }
}
