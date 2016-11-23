package org.broadinstitute.hellbender.utils.bwa;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFlag;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.NativeUtils;

import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.List;

public final class BwaMemIndex implements AutoCloseable {
    private long indexAddress;
    private final List<String> refContigNames;
    private static volatile BwaMemIndex globalInstance;

    public static void createIndexImage( final String refName, final String imgName ) {
        assertNonEmptyReadable(refName+".amb");
        assertNonEmptyReadable(refName+".ann");
        assertNonEmptyReadable(refName+".bwt");
        assertNonEmptyReadable(refName+".pac");
        assertNonEmptyReadable(refName+".sa");
        createIndexImageFile(refName, imgName);
    }

    public static BwaMemIndex getInstance( final String indexImageFile ) {
        BwaMemIndex result = globalInstance;
        if ( result == null ) {
            synchronized(BwaMemIndex.class) {
                result = globalInstance;
                if ( result == null ) globalInstance = result = new BwaMemIndex(indexImageFile);
            }
        }
        return result;
    }

    public static void closeInstance() {
        BwaMemIndex result = globalInstance;
        if ( result != null ) {
            synchronized(BwaMemIndex.class) {
                result = globalInstance;
                if ( result != null ) {
                    globalInstance = null;
                    result.close();
                }
            }
        }
    }

    public static void closeAllDistributedInstances( final JavaSparkContext ctx ) {
        int nJobs = ctx.defaultParallelism();
        final List<Integer> jobList = new ArrayList<>(nJobs);
        for ( int idx = 0; idx != nJobs; ++idx ) jobList.add(idx);
        ctx.parallelize(jobList, nJobs).foreach(idx -> BwaMemIndex.closeInstance());
    }

    static {
        final String libName;
        if (NativeUtils.runningOnLinux()) libName = "/libbwa.Linux.so";
        else if (NativeUtils.runningOnMac()) libName = "/libbwa.Darwin.dylib";
        else libName = null;
        if ( libName == null ) {
            throw new UserException.HardwareFeatureException("We have a JNI binding for bwa-mem only for Linux and Mac.");
        }
        if ( !NativeUtils.loadLibraryFromClasspath(libName) ) {
            throw new UserException.HardwareFeatureException("Misconfiguration: Unable to load bwa-mem native library "+libName);
        }
    }

    public BwaMemIndex( final String indexImageFile ) {
        assertNonEmptyReadable(indexImageFile);
        indexAddress = openIndex(indexImageFile);
        if ( indexAddress == 0L ) {
            throw new GATKException("Unable to open bwa-mem index "+indexImageFile);
        }
        ByteBuffer refContigNamesBuf = getRefContigNames(indexAddress);
        if ( refContigNamesBuf == null ) {
            throw new GATKException("Unable to retrieve reference contig names from bwa-mem index "+indexImageFile);
        }
        refContigNamesBuf.order(ByteOrder.nativeOrder()).position(0).limit(refContigNamesBuf.capacity());
        int nRefContigNames = refContigNamesBuf.getInt();
        refContigNames = new ArrayList<>(nRefContigNames);
        for ( int idx = 0; idx < nRefContigNames; ++idx ) {
            int nameLen = refContigNamesBuf.getInt();
            byte[] nameBytes = new byte[nameLen];
            refContigNamesBuf.get(nameBytes);
            refContigNames.add(new String(nameBytes));
        }
        destroyByteBuffer(refContigNamesBuf);
    }

    @Override
    public synchronized void close() {
        if ( indexAddress != 0 ) destroyIndex(indexAddress);
        indexAddress = 0;
    }

    public BwaMemAligner createAligner() { return new BwaMemAligner(); }

    public List<String> getReferenceContigNames() {
        return refContigNames;
    }

    public final class BwaMemAligner implements Aligner {
        private final ByteBuffer opts;

        public BwaMemAligner() {
            opts = createDefaultOptions();
            opts.order(ByteOrder.nativeOrder()).position(0).limit(opts.capacity());
        }

        public void close() {
            destroyByteBuffer(opts);
        }

        public int getMatchScoreOption() { return opts.getInt(0); }
        public void setMatchScoreOption( final int a ) { opts.putInt(0, a); }
        public int getMismatchPenaltyOption() { return opts.getInt(4); }
        public void setMismatchPenaltyOption( final int b ) { opts.putInt(4, b); }
        public int getDGapOpenPenaltyOption() { return opts.getInt(8); }
        public void setDGapOpenPenaltyOption( final int o_del ) { opts.putInt(8, o_del); }
        public int getDGapExtendPenaltyOption() { return opts.getInt(12); }
        public void setDGapExtendPenaltyOption( final int e_del ) { opts.putInt(12, e_del); }
        public int getIGapOpenPenaltyOption() { return opts.getInt(16); }
        public void setIGapOpenPenaltyOption( final int o_ins ) { opts.putInt(16, o_ins); }
        public int getIGapExtendPenaltyOption() { return opts.getInt(20); }
        public void setIGapExtendPenaltyOption( final int e_ins ) { opts.putInt(20, e_ins); }
        public int getUnpairedPenaltyOption() { return opts.getInt(24); }
        public void setUnpairedPenaltyOption( final int pen_unpaired ) { opts.putInt(24, pen_unpaired); }
        public int getClip5PenaltyOption() { return opts.getInt(28); }
        public void setClip5PenaltyOption( final int pen_clip5 ) { opts.putInt(28, pen_clip5); }
        public int getClip3PenaltyOption() { return opts.getInt(32); }
        public void setClip3PenaltyOption( final int pen_clip3 ) { opts.putInt(32, pen_clip3); }
        public int getBandwidthOption() { return opts.getInt(36); }
        public void setBandwidthOption( final int w ) { opts.putInt(36, w); }
        public int getZDropOption() { return opts.getInt(40); }
        public void setZDropOption( final int zdrop ) { opts.putInt(40, zdrop); }
        public long getMaxMemIntvOption() { return opts.getLong(48); }
        public void setMaxMemIntvOption( final long max_mem_intv ) { opts.putLong(48, max_mem_intv); }
        public int getOutputScoreThresholdOption() { return opts.getInt(56); }
        public void setOutputScoreThresholdOption( final int T ) { opts.putInt(56, T); }

        // flag bits for the flag option
        public static final int MEM_F_PE = 0x2; // this one's particularly important -- the "paired ends" flag
        public static final int MEM_F_NOPAIRING = 0x4;
        public static final int MEM_F_ALL = 0x8;
        public static final int MEM_F_NO_MULTI = 0x10;
        public static final int MEM_F_NO_RESCUE = 0x20;
        public static final int MEM_F_REF_HDR = 0x100;
        public static final int MEM_F_SOFTCLIP = 0x200;
        public static final int MEM_F_SMARTPE = 0x400;
        public static final int MEM_F_PRIMARY5 = 0x800;
        public int getFlagOption() { return opts.getInt(60); }
        public void setFlagOption( final int flag ) { opts.putInt(60, flag); }

        public int getMinSeedLengthOption() { return opts.getInt(64); }
        public void setMinSeedLengthOption( final int min_seed_len ) { opts.putInt(64, min_seed_len); }
        public int getMinChainWeightOption() { return opts.getInt(68); }
        public void setMinChainWeightOption( final int min_chain_weight ) { opts.putInt(68, min_chain_weight); }
        public int getMaxChainExtendOption() { return opts.getInt(72); }
        public void setMaxChainExtendOption( final int max_chain_extend ) { opts.putInt(72, max_chain_extend); }
        public float getSplitFactorOption() { return opts.getFloat(76); }
        public void setSplitFactorOption( final float split_factor ) { opts.putFloat(76, split_factor); }
        public int getSplitWidthOption() { return opts.getInt(80); }
        public void setSplitWidthOption( final int split_width ) { opts.putInt(80, split_width); }
        public int getMaxSeedOccurencesOption() { return opts.getInt(84); }
        public void setMaxSeedOccurencesOption( final int max_occ ) { opts.putInt(84, max_occ); }
        public int getMaxChainGapOption() { return opts.getInt(88); }
        public void setMaxChainGapOption( final int max_chain_gap ) { opts.putInt(88, max_chain_gap); }
        public int getNThreadsOption() { return opts.getInt(92); }
        public void setNThreadsOption( final int n_threads ) { opts.putInt(92, n_threads); }
        public int getChunkSizeOption() { return opts.getInt(96); }
        public void setChunkSizeOption( final int chunk_size ) { opts.putInt(96, chunk_size); }
        public float getMaskLevelOption() { return opts.getFloat(100); }
        public void setMaxLevelOption( final float max_level ) { opts.putFloat(100, max_level); }
        public float getDropRatioOption() { return opts.getFloat(104); }
        public void setDropRatioOption( final float drop_ratio ) { opts.putFloat(104, drop_ratio); }
        public float getXADropRatio() { return opts.getFloat(108); }
        public void setXADropRatio( final float XA_drop_ratio ) { opts.putFloat(108, XA_drop_ratio); }
        public float getMaskLevelRedunOption() { return opts.getFloat(112); }
        public void setMaskLevelRedunOption( final float max_level_redun ) { opts.putFloat(112, max_level_redun); }
        public float getMapQCoefLenOption() { return opts.getFloat(116); }
        public void setMapQCoefLenOption( final float mapQ_coef_len ) { opts.putFloat(116, mapQ_coef_len); }
        public int getMapQCoefFacOption() { return opts.getInt(120); }
        public void setMapQCoefFacOption( final int mapQ_coef_fac ) { opts.putInt(120, mapQ_coef_fac); }
        public int getMaxInsOption() { return opts.getInt(124); }
        public void setMaxInsOption( final int max_ins ) { opts.putInt(124, max_ins); }
        public int getMaxMateSWOption() { return opts.getInt(128); }
        public void setMaxMateSWOption( final int max_matesw ) { opts.putInt(128, max_matesw); }
        public int getMaxXAHitsOption() { return opts.getInt(132); }
        public void setMaxXAHitsOption( final int max_XA_hits ) { opts.putInt(132, max_XA_hits); }
        public int getMaxXAHitsAltOption() { return opts.getInt(136); }
        public void setMaxXAHitsAltOption( final int max_XA_hits_alt ) { opts.putInt(136, max_XA_hits_alt); }
        public byte[] getScoringMatrixOption() {
            final byte[] result = new byte[25];
            opts.position(140);
            opts.get(result);
            return result; }
        public void setScoringMatrixOption( final byte[] mat ) {
            opts.position(140);
            opts.put(mat);
        }

        public void setIntraCtgOptions() {
            setDGapOpenPenaltyOption(16);
            setIGapOpenPenaltyOption(16);
            setMismatchPenaltyOption(9);
            setClip5PenaltyOption(5);
            setClip3PenaltyOption(5);
        }

        public List<String> getReferenceContigNames() {
            return refContigNames;
        }

        public List<List<Alignment>> alignSeqs( List<byte[]> sequences ) {
            if ( indexAddress == 0L ) {
                throw new GATKException("No bwa-mem index is open.");
            }
            final int nContigs = sequences.size();
            // 4 bytes for the initial contig count, a null byte at the end of each sequence, and all the sequence bytes
            final int capacity = 4 + nContigs + sequences.stream().mapToInt(seq -> seq.length).sum();
            final ByteBuffer contigBuf = ByteBuffer.allocateDirect(capacity);
            contigBuf.order(ByteOrder.nativeOrder());
            contigBuf.putInt(nContigs);
            sequences.forEach(seq -> contigBuf.put(seq).put((byte)0));
            contigBuf.flip();
            ByteBuffer alignsBuf = createAlignments(contigBuf, indexAddress, opts);
            if ( alignsBuf == null ) {
                throw new GATKException("Unable to get alignments from bwa-mem. We don't know why.");
            }
            alignsBuf.order(ByteOrder.nativeOrder()).position(0).limit(alignsBuf.capacity());
            final List<List<Alignment>> allAlignments = new ArrayList<>(nContigs);
            for ( int contigId = 0; contigId != nContigs; ++contigId ) {
                final int contigLength = sequences.get(contigId).length;
                int nAligns = alignsBuf.getInt();
                final List<Alignment> alignments = new ArrayList<>(nAligns);
                while ( nAligns-- > 0 ) {
                    final int flag_mapQ = alignsBuf.getInt();
                    final int flags = flag_mapQ >>> 16;
                    final int mapQual = flag_mapQ & 0xff;
                    final int refId;
                    final int refStart;
                    final int refEnd;
                    final int tigStart;
                    final int tigEnd;
                    final int nMismatches;
                    final int alignerScore;
                    final int suboptimalScore;
                    final Cigar cigar = new Cigar();
                    final String mdTag;
                    final String xaTag;

                    // if unmapped
                    if ( (flags & SAMFlag.READ_UNMAPPED.intValue()) != 0 ) {
                        refId = -1;
                        refStart = -1;
                        refEnd = -1;
                        tigStart = -1;
                        tigEnd = -1;
                        nMismatches = 0;
                        alignerScore = 0;
                        suboptimalScore = 0;
                        mdTag = null;
                        xaTag = null;
                    }
                    else { // mapped
                        refId = alignsBuf.getInt();
                        refStart = alignsBuf.getInt();
                        nMismatches = alignsBuf.getInt();
                        alignerScore = alignsBuf.getInt();
                        suboptimalScore = alignsBuf.getInt();
                        int nCigarOps = alignsBuf.getInt();
                        while ( nCigarOps-- > 0 ) {
                            final int lenOp = alignsBuf.getInt();
                            cigar.add(new CigarElement(lenOp >>> 4, CigarOperator.binaryToEnum(lenOp & 0x0f)));
                        }
                        mdTag = getTag(alignsBuf);
                        xaTag = getTag(alignsBuf);
                        refEnd = refStart + cigar.getReferenceLength();
                        final CigarElement ele0 = cigar.getFirstCigarElement();
                        tigStart = ele0.getOperator() == CigarOperator.S ? ele0.getLength() : 0;
                        final CigarElement eleN = cigar.getLastCigarElement();
                        tigEnd = contigLength - (eleN.getOperator() == CigarOperator.S ? eleN.getLength() : 0);
                    }

                    final int mateRefId;
                    final int mateStartPos;
                    final int templateLen;

                    // if unpaired, or mate unmapped
                    if ( (flags & SAMFlag.READ_PAIRED.intValue()) == 0 ||
                            (flags & SAMFlag.MATE_UNMAPPED.intValue()) != 0 ) {
                        mateRefId = -1;
                        mateStartPos = -1;
                        templateLen = 0;
                    } else { // has mapped mate
                        mateRefId = alignsBuf.getInt();
                        mateStartPos = alignsBuf.getInt();
                        templateLen = alignsBuf.getInt();
                    }
                    alignments.add(new Alignment(flags, refId, refStart, refEnd, tigStart, tigEnd, mapQual,
                                                 nMismatches, alignerScore, suboptimalScore, cigar, mdTag, xaTag,
                                                 mateRefId, mateStartPos, templateLen));
                }
                allAlignments.add(alignments);
            }
            destroyByteBuffer(alignsBuf);
            return allAlignments;
        }

        private String getTag( final ByteBuffer buffer ) {
            int tagLen = buffer.getInt();
            if ( tagLen == 0 ) return null;
            byte[] tagBytes = new byte[(tagLen+3)&~3];
            buffer.get(tagBytes);
            return new String(tagBytes, 0, tagLen);
        }
    }

    private static void assertNonEmptyReadable( final String fileName ) {
        if ( !nonEmptyReadableFile(fileName) )
            throw new GATKException("Missing bwa index file: "+fileName);
    }

    private static boolean nonEmptyReadableFile( final String fileName ) {
        try ( final FileInputStream is = new FileInputStream(fileName) ) {
            return is.read() != -1;
        } catch ( final IOException ioe ) {
            return false;
        }
    }

    private static native boolean createIndexImageFile( String referenceName, String imageName );
    private static native long openIndex( String indexImageFile );
    private static native int destroyIndex( long indexAddress );
    private static native ByteBuffer createDefaultOptions();
    private static native ByteBuffer getRefContigNames( long indexAddress );
    private static native ByteBuffer createAlignments( ByteBuffer seqs, long indexAddress, ByteBuffer opts );
    private static native void destroyByteBuffer( ByteBuffer alignments );
    public static native String getVersion();
}
