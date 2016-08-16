/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.sam.markduplicates;

import htsjdk.samtools.DuplicateSet;
import htsjdk.samtools.DuplicateSetIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;

import java.util.*;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.counting;

public class UmiAwareDuplicateSetIterator extends UmiGraph implements CloseableIterator<DuplicateSet> {
    private final DuplicateSetIterator wrappedIterator;
    private Iterator<DuplicateSet> nextSetsIterator;
    private final int editDistanceToJoin;
    private final boolean addInferredUmi;
    private final String umiTag;
    private final String inferredUmiTag;

    public UmiAwareDuplicateSetIterator(final DuplicateSetIterator wrappedIterator, final int editDistanceToJoin, final boolean addInferredUmi, final String umiTag, final String inferredUmiTag) {
        this.wrappedIterator = wrappedIterator;
        this.editDistanceToJoin = editDistanceToJoin;
        this.addInferredUmi = addInferredUmi;
        this.umiTag = umiTag;
        this.inferredUmiTag = inferredUmiTag;
        nextSetsIterator = Collections.emptyIterator();
    }

    @Override
    public void close() {
        wrappedIterator.close();
    }

    @Override
    public boolean hasNext() {
        return nextSetsIterator.hasNext() || wrappedIterator.hasNext();
    }

    @Override
    public DuplicateSet next() {
        if(!nextSetsIterator.hasNext())
            process(wrappedIterator.next());

        return nextSetsIterator.next();
    }

    // Takes a duplicate set and breaks it up into possible smaller sets according to the UMI,
    // and updates nextSetsIterator to be an iterator on that set of DuplicateSets.
    private void process(final DuplicateSet set) {

        final List<SAMRecord> records = set.getRecords();

        // If any records are missing the UMI_TAG proceed as if there were no UMIs
        // and return nextSetsIterator without breaking it up into smaller sets.
        for(SAMRecord rec : records) {
            if(rec.getStringAttribute(umiTag) == null) {
                nextSetsIterator = Collections.singleton(set).iterator();
                return;
            }
        }

        // Count all the uniquely observed UMI sequences
        Map<String, Long> umiCounts = records.stream()
                        .collect(Collectors.groupingBy(p -> p.getStringAttribute(umiTag),
                                counting()));

        UmiGraph umiGraph = new UmiGraph(umiCounts, editDistanceToJoin);

        List<DuplicateSet> duplicateSetList = umiGraph.joinUmisIntoDuplicateSets(records, umiTag, inferredUmiTag);

        nextSetsIterator = duplicateSetList.iterator();
    }


}

