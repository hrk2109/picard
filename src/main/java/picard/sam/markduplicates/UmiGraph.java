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
import htsjdk.samtools.SAMRecord;
import picard.PicardException;

import java.util.*;

public class UmiGraph {
    private Map<String, Long> umiCounts;
    private int[] duplicateSet;    // id[i] = parent of i
    private Long[] occupancy;
    private String[] umi;
    private int editDistanceToJoin;
    private int nUmis;
    private int nDuplicateSets;

    UmiGraph() {

    }

    UmiGraph(final Map<String, Long> umiCounts, final int editDistanceToJoin) {
        // At first every umi is it's own duplicate set
        nUmis = umiCounts.size();
        nDuplicateSets = nUmis;

        occupancy = new Long[nUmis];
        umi = new String[nUmis];
        this.umiCounts = umiCounts;
        this.editDistanceToJoin = editDistanceToJoin;

        duplicateSet = new int[nUmis];
        for (int i = 0; i < nUmis; i++) {
            duplicateSet[i] = i;
        }

        int i = 0;
        for (Map.Entry<String, Long> entry : umiCounts.entrySet()) {
            umi[i] = entry.getKey();
            occupancy[i] = entry.getValue();
            i++;
        }
    }

    public int find(int p) {
        int root = p;
        while (root != duplicateSet[root])
            root = duplicateSet[root];
        while (p != root) {
            int newp = duplicateSet[p];
            duplicateSet[p] = root;
            p = newp;
        }
        return root;
    }

    public void union(int p, int q) {
        int rootP = find(p);
        int rootQ = find(q);
        if (rootP == rootQ) return;
        duplicateSet[rootP] = rootQ;
        nDuplicateSets--;
    }

    List<DuplicateSet> joinUmisIntoDuplicateSets(List<SAMRecord> records, String umiTag, String inferredUmiTag) {
        // Iterate over diagonal
        for (int i = 0; i < nUmis; i++) {
            for (int j = i + 1; j < nUmis; j++) {
                if (getEditDistance(umi[i], umi[j]) <= editDistanceToJoin) {
                    union(i, j);
                }
            }
        }

        for (int i = 0; i < nUmis; i++) {
            duplicateSet[i] = find(i);
        }

        Map<Integer, List<SAMRecord>> duplicateSets = new HashMap<>();

        // Assign UMIs to duplicateSets
        Map<String, Integer> duplicateSetsFromUmis = getDuplicateSetsFromUmis();
        for(SAMRecord rec : records) {

            String umi = rec.getStringAttribute(umiTag);
            Integer duplicateSetIndex = duplicateSetsFromUmis.get(umi);

            if(duplicateSets.containsKey(duplicateSetIndex)) {
                duplicateSets.get(duplicateSetIndex).add(rec);
            }
            else {
                List<SAMRecord> n = new ArrayList();
                n.add(rec);
                duplicateSets.put(duplicateSetIndex, n);
            }
        }


        List<DuplicateSet> duplicateSetList = new ArrayList();
        for (Map.Entry<Integer, List<SAMRecord>> entry : duplicateSets.entrySet()) {
            DuplicateSet ds = new DuplicateSet();
            List<SAMRecord> recordList = entry.getValue();

            // Add records to the DuplicateSet
            for(SAMRecord rec : recordList) {
                ds.add(rec);
            }

            // For a particular duplicate set, identify the most common umi
            // and use this an an iferred umi
            Long maxCount = new Long(0);
            String inferredUmi = null;
            for(SAMRecord rec : recordList) {
                String umi = rec.getStringAttribute(umiTag);
                //inferredUmi = umi;
                if(umiCounts.get(umi) > maxCount) {
                   maxCount = umiCounts.get(umi);
                   inferredUmi = umi;
                }
            }

            // Set the records to contain the inferred Umi
            for(SAMRecord rec : recordList) {
                rec.setAttribute(inferredUmiTag, inferredUmi);
            }

            duplicateSetList.add(ds);
        }

        return duplicateSetList;
    }

    private int getEditDistance(final String s1, final String s2) {
        if(s1 == null || s2 == null) {
            throw new PicardException("Attempt to compare two incomparable UMIs.  At least one of the UMIs was null.");
        }
        if(s1.length() != s2.length()) {
            throw new PicardException("Barcode " + s1 + " and " + s2 + " do not have matching lengths.");
        }
        int count = 0;
        for(int i = 0;i < s1.length();i++) {
            if(s1.charAt(i) != s2.charAt(i)) {
                count++;
            }
        }
        return count;
    }

    Map<Integer, List<String>> getUmisFromDuplicateSet() {
        Map<Integer, List<String>> umisByDuplicateSet = new HashMap<>();

        for(int i = 0;i < duplicateSet.length;i++) {
            if(!umisByDuplicateSet.containsKey(duplicateSet[i])) {
                umisByDuplicateSet.put(duplicateSet[i], new ArrayList<>());
            }
            umisByDuplicateSet.get(duplicateSet[i]).add(umi[i]);
        }

        return umisByDuplicateSet;
    }

    Map<String, Integer> getDuplicateSetsFromUmis() {
        Map<String, Integer> duplicateSetsFromUmis = new HashMap<>();
        for(int i = 0;i < duplicateSet.length;i++) {
            duplicateSetsFromUmis.put(umi[i], duplicateSet[i]);
        }
        return duplicateSetsFromUmis;
    }

    public String getInferredUmi(String umi) {

        Map<Integer, List<String>> umisFromDuplicateSet = getUmisFromDuplicateSet();
        Map<String, Integer> duplicateSetsFromUmis = getDuplicateSetsFromUmis();

        List<String> umisInDuplicateSet = umisFromDuplicateSet.get(duplicateSetsFromUmis.get(umi));

        // Find the most commonly occuring umi in the duplicate set


        for(String umiInDuplicateSet : umisInDuplicateSet) {
            umiCounts.get(umiInDuplicateSet);
        }

        return umi;
    }
}
