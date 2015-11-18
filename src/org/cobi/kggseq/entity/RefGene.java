/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author MX Li
 */
public class RefGene implements Cloneable, Serializable {

    protected String symb;
    protected int entrezID = -1;
    protected String fullName;
    protected int oMIMID;
    protected int length;
    protected Map<String, RefRefmRNA> mRNAList;
    protected String chromosome;
    protected int startPos;
    protected int endPos;
    protected int order;
    protected String description;
    protected List<SeqSegment> mergedSegments;
    protected int totalMergedLen = 0;
    protected char strand;

    public int getTotalMergedLen() {
        return totalMergedLen;
    }

    public void setTotalMergedLen(int totalMergedLen) {
        this.totalMergedLen = totalMergedLen;
    }

    public char getStrand() {
        return strand;
    }

    public void setStrand(char strand) {
        this.strand = strand;
    }

    public List<SeqSegment> getMergedSegments() {
        return mergedSegments;
    }

    public void setExtractedSegments(List<SeqSegment> extractedSegments) {
        this.mergedSegments = extractedSegments;
    }

    //as a gene may have various transcripts, they may have different exons
    public void combineExons(int exonExtendLen, boolean onlyCoding) throws Exception {
        if (mergedSegments == null) {
            mergedSegments = new ArrayList<SeqSegment>();
        }
        boolean isForward = true;
        if (strand == '-') {
            isForward = false;
        }
        // System.out.println(this.getSymb());
   
        for (Map.Entry<String, RefRefmRNA> m2 : mRNAList.entrySet()) {
            RefRefmRNA mrna = m2.getValue();
            List<RefExon> exonList = mrna.getExons();
            mergeExons2Segs(exonList, exonExtendLen, mergedSegments, onlyCoding, isForward);
        }
    }

    //as a gene may have various transcripts, they may have different exons
    public void combineIntron(int intronNum) throws Exception {
        if (intronNum <= 0) {
            return;
        }
        Map<String, RefRefmRNA> mRNAMap = this.getMRNAList();

        if (mergedSegments == null) {
            mergedSegments = new ArrayList<SeqSegment>();
        }

        for (Map.Entry<String, RefRefmRNA> m2 : mRNAMap.entrySet()) {
            RefRefmRNA mrna = m2.getValue();
            List<RefIntron> intronList = mrna.getIntrons();
            mergeIntrons2Segs(intronList, intronNum, mergedSegments);
        }
    }

    public void mergeIntrons2Segs(final List<RefIntron> origIntronList, int maxNum, List<SeqSegment> destSegList) throws Exception {
        int intronIndex = 0;
        int segIndex = 0;
        List<SeqSegment> tmpSegList = new ArrayList<SeqSegment>();
        tmpSegList.addAll(destSegList);
        destSegList.clear();
        int nIntron = origIntronList.size();
        int nSeg = tmpSegList.size();
        int intronStart, intronEnd, segStart, segEnd;
        //note the origIntronList and destSegList must be sorted
        while ((intronIndex < maxNum) && (intronIndex < nIntron) && (segIndex < nSeg)) {
            intronStart = origIntronList.get(intronIndex).getStart();
            intronEnd = origIntronList.get(intronIndex).getEnd();
            segStart = tmpSegList.get(segIndex).getStart();
            segEnd = tmpSegList.get(segIndex).getEnd();
            //six possible relative positions
            if (segEnd < intronStart) {
                destSegList.add(tmpSegList.get(segIndex));
                segIndex++;
            } else if (segEnd >= intronStart && segEnd <= intronEnd && segStart <= intronStart) {
                //merged the two and move intronIndex
                tmpSegList.get(segIndex).setEnd(intronEnd);
                tmpSegList.get(segIndex).setDescription(tmpSegList.get(segIndex).getDescription() + "+Intron");
                intronIndex++;
            } else if (segStart <= intronStart && segEnd >= intronEnd) {
                //overlapped, move seg index
                intronIndex++;
                tmpSegList.get(segIndex).setDescription(tmpSegList.get(segIndex).getDescription() + "+Intron");
            } else if (intronEnd < segStart) {
                //create a new segment for this Exon and move intronIndex
                SeqSegment seg = new SeqSegment();
                seg.setStart(intronStart);
                seg.setEnd(intronEnd);
                seg.setDescription("Intron");
                destSegList.add(seg);
                intronIndex++;
            } else if (segStart <= intronEnd && segStart >= intronStart && segEnd >= intronEnd) {
                //merged the two and move seg index
                tmpSegList.get(segIndex).setStart(intronStart);
                tmpSegList.get(segIndex).setDescription(tmpSegList.get(segIndex).getDescription() + "+Intron");
                intronIndex++;
            } else if (segStart >= intronStart && segEnd <= intronEnd) {
                //overlapped, move seg index
                tmpSegList.get(segIndex).setStart(intronStart);
                tmpSegList.get(segIndex).setEnd(intronEnd);
                tmpSegList.get(segIndex).setDescription(tmpSegList.get(segIndex).getDescription() + "+Intron");
                intronIndex++;
            }
        }

        while (intronIndex < maxNum && intronIndex < nIntron) {
            intronStart = origIntronList.get(intronIndex).getStart();
            intronEnd = origIntronList.get(intronIndex).getEnd();
            SeqSegment seg = new SeqSegment();
            seg.setStart(intronStart);
            seg.setEnd(intronEnd);
            seg.setDescription("Intron");
            destSegList.add(seg);
            intronIndex++;
        }
        while (segIndex < nSeg) {
            destSegList.add(tmpSegList.get(segIndex));
            segIndex++;
        }

    }

    public void mergeExons2Segs(final List<RefExon> origExonList, int exonExtendLen, List<SeqSegment> destSegList, boolean onlyCoding, boolean isForward) throws Exception {
        int exonIndex = 0;
        int segIndex = 0;
        List<SeqSegment> tmpSegList = new ArrayList<SeqSegment>();
        tmpSegList.addAll(destSegList);
        destSegList.clear();
        int nExon = origExonList.size();
        int nSeg = tmpSegList.size();
        int exonStart, exonEnd, segStart, segEnd;
        //note the origIntronList and destSegList must be sorted
        while ((exonIndex < nExon) && (segIndex < nSeg)) {
            if (onlyCoding) {
                if (origExonList.get(exonIndex).getCodingLength() <= 0) {
                    exonIndex++;
                    continue;
                }
                exonStart = origExonList.get(exonIndex).getCodingStart() - exonExtendLen;
                exonEnd = origExonList.get(exonIndex).getCodingEnd() + exonExtendLen;
            } else {
                exonStart = origExonList.get(exonIndex).getStart() - exonExtendLen;
                exonEnd = origExonList.get(exonIndex).getEnd() + exonExtendLen;
            }

            segStart = tmpSegList.get(segIndex).getStart();
            segEnd = tmpSegList.get(segIndex).getEnd();
            //start with 3' comparision
            //seven possible relative positions
            //assume all start are smaller than the end regardless of the strand
            if (isForward) {
                if (segStart > exonEnd) {
                    SeqSegment seg = new SeqSegment();
                    seg.setStart(exonStart);
                    seg.setEnd(exonEnd);
                    seg.setDescription("Exon");
                    destSegList.add(seg);
                    exonIndex++;
                } else {
                    if (segStart >= exonStart && segEnd >= exonEnd) {
                        tmpSegList.get(segIndex).setStart(exonStart);
                        exonIndex++;
                        tmpSegList.get(segIndex).setDescription(tmpSegList.get(segIndex).getDescription() + "+Exon");
                    } else if (segStart < exonStart && segEnd >= exonEnd) {
                        exonIndex++;
                        tmpSegList.get(segIndex).setDescription(tmpSegList.get(segIndex).getDescription() + "+Exon");
                    } else if (segStart >= exonStart && segEnd < exonEnd) {
                        tmpSegList.get(segIndex).setStart(exonStart);
                        tmpSegList.get(segIndex).setEnd(exonEnd);
                        exonIndex++;
                        tmpSegList.get(segIndex).setDescription(tmpSegList.get(segIndex).getDescription() + "+Exon");
                    } else if (segEnd < exonStart) {
                        destSegList.add(tmpSegList.get(segIndex));
                        segIndex++;
                    } else if (segStart < exonStart && segEnd <= exonEnd) {
                        exonIndex++;
                        tmpSegList.get(segIndex).setDescription(tmpSegList.get(segIndex).getDescription() + "+Exon");
                        tmpSegList.get(segIndex).setEnd(exonEnd);
                    }
                }
            } else {
                if (segEnd < exonStart) {
                    //a new region
                    SeqSegment seg = new SeqSegment();
                    seg.setStart(exonStart);
                    seg.setEnd(exonEnd);
                    seg.setDescription("Exon");
                    destSegList.add(seg);
                    exonIndex++;
                } else {
                    if (segStart <= exonStart && segEnd <= exonEnd) {
                        tmpSegList.get(segIndex).setEnd(exonEnd);
                        tmpSegList.get(segIndex).setDescription(tmpSegList.get(segIndex).getDescription() + "+Exon");
                        exonIndex++;
                    } else if (segStart <= exonStart && segEnd > exonEnd) {
                        tmpSegList.get(segIndex).setDescription(tmpSegList.get(segIndex).getDescription() + "+Exon");
                        exonIndex++;
                    } else if (segStart > exonStart && segEnd <= exonEnd) {
                        tmpSegList.get(segIndex).setStart(exonStart);
                        tmpSegList.get(segIndex).setEnd(exonEnd);
                        exonIndex++;
                        tmpSegList.get(segIndex).setDescription(tmpSegList.get(segIndex).getDescription() + "+Exon");
                    } else if (segStart > exonEnd) {
                        destSegList.add(tmpSegList.get(segIndex));
                        segIndex++;
                    } else if (segStart > exonStart && segEnd > exonEnd) {
                        tmpSegList.get(segIndex).setStart(exonStart);
                        tmpSegList.get(segIndex).setDescription(tmpSegList.get(segIndex).getDescription() + "+Exon");
                        exonIndex++;
                    }
                }
            }
        }

        while (exonIndex < nExon) {
            if (onlyCoding) {
                if (origExonList.get(exonIndex).getCodingLength() <= 0) {
                    exonIndex++;
                    continue;
                }
                exonStart = origExonList.get(exonIndex).getCodingStart() - exonExtendLen;
                exonEnd = origExonList.get(exonIndex).getCodingEnd() + exonExtendLen;
            } else {
                exonStart = origExonList.get(exonIndex).getStart() - exonExtendLen;
                exonEnd = origExonList.get(exonIndex).getEnd() + exonExtendLen;
            }

            SeqSegment seg = new SeqSegment();
            seg.setStart(exonStart);
            seg.setEnd(exonEnd);
            seg.setDescription("Exon");
            destSegList.add(seg);
            exonIndex++;
        }
        while (segIndex < nSeg) {
            destSegList.add(tmpSegList.get(segIndex));
            segIndex++;
        }
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        RefGene o = (RefGene) super.clone();
        o.mergedSegments = new ArrayList<SeqSegment>();
        return o;
    }

    public void writeObject(ObjectOutputStream outputStream) throws IOException {
        outputStream.defaultWriteObject();
        // outputStream.writeObject(mRNAList);
    }

    private void readObject(ObjectInputStream inputStream) throws IOException, ClassNotFoundException {
        inputStream.defaultReadObject();
        //mRNAList = (Map<String, mRNA>) inputStream.readObject();
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public int getOrder() {
        return order;
    }

    public void setOrder(int order) {
        this.order = order;
    }

    public int getEndPos() {
        return endPos;
    }

    public void setEndPos(int endPos) {
        this.endPos = endPos;
    }

    public int getStartPos() {
        return startPos;
    }

    public void setStartPos(int startPos) {
        this.startPos = startPos;
    }

    public String getChromosome() {
        return chromosome;
    }

    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    public int getLength() {
        return length;
    }

    public void setLength(int length) {
        this.length = length;
    }

    public Map<String, RefRefmRNA> getMRNAList() {
        return mRNAList;
    }

    public void setMRNAList(Map<String, RefRefmRNA> mRNAList) {
        this.mRNAList = mRNAList;
    }

    public RefGene(String symb, int entrezID) {
        this.symb = symb;
        this.entrezID = entrezID;
    }

    public void addmRNA(RefRefmRNA mrna) {
        if (mRNAList == null) {
            mRNAList = new HashMap<String, RefRefmRNA>();
        }
        mRNAList.put(mrna.getRefID(), mrna);
    }

    public RefGene() {
    }

    public RefGene(String symb) {
        this.symb = symb;
    }

    /**
     * Get the value of oMIMID
     *
     * @return the value of oMIMID
     */
    public int getOMIMID() {
        return oMIMID;
    }

    /**
     * Set the value of oMIMID
     *
     * @param oMIMID new value of oMIMID
     */
    public void setOMIMID(int oMIMID) {
        this.oMIMID = oMIMID;
    }

    /**
     * Get the value of fullName
     *
     * @return the value of fullName
     */
    public String getFullName() {
        return fullName;
    }

    /**
     * Set the value of fullName
     *
     * @param fullName new value of fullName
     */
    public void setFullName(String fullName) {
        this.fullName = fullName;
    }

    /**
     * Get the value of entrezID
     *
     * @return the value of entrezID
     */
    public int getEntrezID() {
        return entrezID;
    }

    /**
     * Set the value of entrezID
     *
     * @param entrezID new value of entrezID
     */
    public void setEntrezID(int entrezID) {
        this.entrezID = entrezID;
    }

    /**
     * Get the value of symb
     *
     * @return the value of symb
     */
    public String getSymb() {
        return symb;
    }

    /**
     * Set the value of symb
     *
     * @param symb new value of symb
     */
    public void setSymb(String symb) {
        this.symb = symb;
    }
}
