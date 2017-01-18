package org.broadinstitute.hellbender.tools.spark.bwa;

import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import static org.testng.Assert.assertEquals;

public class BwaSparkEngineUnitTest extends BaseTest {

    @DataProvider(name="readPairsAndPartitions")
    public Object[][] readPairsAndPartitions() {
        return new Object[][] {
                { 1, 1, new int[] {2} },
                { 2, 2, new int[] {2, 2} },
                { 3, 2, new int[] {4, 2} },
                { 3, 3, new int[] {2, 2, 2} },
                { 6, 2, new int[] {6, 6} },
                { 6, 3, new int[] {4, 4, 4} },
                { 6, 4, new int[] {4, 2, 4, 2} },
        };
    }

    @Test(dataProvider = "readPairsAndPartitions")
    public void testPutPairsInSamePartition(int numPairs, int numPartitions, int[] expectedReadsPerPartition) throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        JavaRDD<GATKRead> reads = createReads(ctx, numPairs, numPartitions);
        JavaRDD<GATKRead> pairedReads = BwaSparkEngine.putPairsInSamePartition(ctx, reads);
        List<List<GATKRead>> partitions = pairedReads.mapPartitions((FlatMapFunction<Iterator<GATKRead>, List<GATKRead>>) it ->
                Iterators.singletonIterator(Lists.newArrayList(it))).collect();
        assertEquals(partitions.size(), numPartitions);
        for (int i = 0; i < numPartitions; i++) {
            assertEquals(partitions.get(i).size(), expectedReadsPerPartition[i]);
        }
        assertEquals(Arrays.stream(expectedReadsPerPartition).sum(), numPairs * 2);
    }

    private JavaRDD<GATKRead> createReads(JavaSparkContext ctx, int numPairs, int numPartitions) {
        SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final int readSize = 151;
        final int fragmentLen = 400;
        final String templateName = "readpair";
        int leftStart = 10000;
        List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < numPairs;i++) {
            leftStart += readSize * 2;
            int rightStart = leftStart + fragmentLen - readSize;
            reads.addAll(ArtificialReadUtils.createPair(header, templateName + i, readSize, leftStart, rightStart, true, false));
        }
        return ctx.parallelize(reads, numPartitions);
    }
}
