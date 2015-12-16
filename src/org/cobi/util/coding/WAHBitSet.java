/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.coding;

import cern.colt.list.IntArrayList;

import java.util.Arrays;
import java.util.Iterator;
import java.util.BitSet;
import org.apache.log4j.Logger;
import java.io.Serializable;

/**
 * WAHBitSet implements the Word-Aligned Hybrid compressed BitSet.
 *
 * Using WAH as the compression scheme the class achieves both compressaion and
 * performance on bitset operations.
 */
public class WAHBitSet {

    private static Logger logger = Logger.getLogger(WAHBitSet.class.getName());

    private static final long serialVersionUID = -7210214822424094012L;

    // constances defined in bitvector
    private static final int MAXBITS = 8 * 4 - 1;
    private static final int SECONDBIT = MAXBITS - 1;
    private static final int ALLONES = Integer.MAX_VALUE;
    private static final int MAXCNT = ((1 << SECONDBIT) - 1);
    private static final int FILLBIT = (1 << SECONDBIT);
    private static final int HEADER0 = (2 << SECONDBIT);
    private static final int HEADER1 = (3 << SECONDBIT);

    private static final boolean RUN_UNTESTED_CODE = false;
    private static final int DEFAULT_INITIAL_SIZE = 4;

    private int nset;
    private int nbits;
    IntArrayList vec = new IntArrayList(DEFAULT_INITIAL_SIZE);
    ActiveWord active = new ActiveWord();

    /**
     * Create an empty bitset.
     */
    public WAHBitSet() {
    }

    /**
     * Initialize the compressed bitset withn an uncompressed bitset.
     *
     * @param set the uncompressed bitset.
     */
    public WAHBitSet(BitSet set) {
        int last = 0, next, size = 0;
        while ((next = set.nextSetBit(last)) != -1) {
            set(next);
            last = next + 1;
            size++;
        }

        if (size != set.cardinality() || size != this.cardinality()) {
            throw new AssertionError("The sizes were not the same after and before!");
        }
    }

    /**
     * Sets the index <tt>i</tt> to one. expands the bitset if required.
     *
     * Note: currently we are going to support only increasing the indexes.
     *
     * @param i the index to set as 1.
     */
    public void set(int i) {
        setBit(i, 1);
    }

    /**
     * Checks if the bit is set in the compressed bitset.
     *
     * This operation is considerably slower than the uncompressed version. Use
     * it with care.
     *
     * @param i the index where bit is checked.
     * @return true if the index is set, false otherwise.
     */
    public boolean get(int i) {
        if (i > numBits()) { // no-way
            return false;
        } else { // we have to do hard way

            for (int j = 0; j < vec.size() && i >= 0; j++) {
                int v = vec.getQuick(j);

                // for a fill, we know easily.
                if (isAFill(v)) {
                    int size = (v & MAXCNT) * MAXBITS;
                    if (i < size) {
                        return isOneFill(v);
                    }
                    i -= size;
                } else // plain value.
                 if (i < MAXBITS) { // do we need to check?
                        return (1 << (MAXBITS - i - 1) & v) != 0;
                    } else { // yet to get there.
                        i -= MAXBITS;
                    }
            }

            // we need to check the active word.
            if (i >= 0) {
                return (moreThanInUnsigned(active.val << (MAXBITS + 1 - active.nbits + i), ALLONES));
            }
        }

        return false;
    }

    /**
     * Returns the number of 1 bits in the bitset.
     *
     * @return the number of 1 bits in the bitset.
     */
    public int cardinality() {
        // check the cardinality, again if it has been reset.
        if (nset == 0 && !vec.isEmpty()) {
            doCount();
        }

        // the sizes in the vector and the active word.
        return nset + Integer.bitCount(active.val);
    }

    /**
     * Returns an iterator for the set bits in the bitset.
     *
     * @return an iterator for the set bits in the bitset
     */
    public Iterator iterator() {
        return new WAHIterator();
    }

    /**
     * Returns an index set for the set bits in the bitset. IndexSet is an
     * easier way to get the set bits in the compressed bitset. IndexSet is
     * lighter weight interface than iterator, and should be used when the
     * performance is important.
     *
     * @return an index set for the set bits in the bitset
     */
    public IndexSet getIndexSet() {
        return new IndexSet();
    }

    /**
     * Returns the amount of memory used by the compressed bit set
     *
     * @return the amount of memory used by the compressed bit set
     */
    public long memSize() {
        return vec.size() + 2;
    }

    /**
     * Returns a new WAH compressed bitset after anding the current bitset with
     * the <i>other</i> bitset.
     *
     * @param other the bitset to and with
     * @return The resulting bitset
     */
    public WAHBitSet and(WAHBitSet other) {
        WAHBitSet ret = new WAHBitSet();

        // ensure that they have the same bit length.
        if (this.numBits() < other.numBits()) {
            this.setBit(other.numBits() - 1, 0);
        } else if (this.numBits() > other.numBits()) {
            other.setBit(numBits() - 1, 0);
        }

        // if there is something in the vector.
        if (vec.size() > 0) {
            // create new run objects and decode them.
            run xrun = new run(vec), yrun = new run(other.vec);
            xrun.decode();
            yrun.decode();
            do {
                // if you finished a run, then get the next one.
                if (xrun.nWords == 0) {
                    xrun.inc();
                    xrun.decode();
                }

                if (yrun.nWords == 0) {
                    yrun.inc();
                    yrun.decode();
                }

                if (xrun.isFill()) {
                    if (yrun.isFill()) {
                        // both are fills... this is the best.
                        int nWords = Math.min(xrun.nWords, yrun.nWords);
                        ret.appendFill(nWords, xrun.fillWord & yrun.fillWord);
                        xrun.nWords -= nWords;
                        yrun.nWords -= nWords;
                    } else {
                        // just cut through the other run
                        chewUpRun(xrun, ret, yrun);
                    }
                } else if (yrun.isFill()) {
                    // again do the same, with different order.
                    chewUpRun(yrun, ret, xrun);
                } else {
                    // both are literals, so get the new literal and
                    // append it to the return value.
                    ret.active.val = xrun.get() & yrun.get();
                    ret.appendLiteral();
                    yrun.nWords = 0;
                    xrun.nWords = 0;
                }
                // till they are not at the end.
            } while (!(xrun.end() && yrun.end()));
        }

        // set the active word.
        ret.active.val = this.active.val & other.active.val;
        ret.active.nbits = this.active.nbits;

        // ensure that the counts are set properly.
        ret.doCount();

        return ret;
    }

    /**
     * Returns a new WAH compressed bitset after oring the current bitset with
     * the <i>other</i> bitset. This function is not as optimized as <i>and</i>
     *
     * @param other the bitset to or with
     * @return The resulting bitset
     */
    public WAHBitSet or(WAHBitSet other) {
        return genericOp(OpType.Or, other);
    }

    private WAHBitSet andNot(WAHBitSet other) {
        return genericOp(OpType.AndNot, other);
    }

    /**
     * This is an optimization over the and function. This does not create new
     * bitset. This just counts the number of 1 bits common between the two
     * bitsets.
     *
     * @param other the bitset to and with.
     * @return the number of 1s common between two bitsets.
     */
    public int andSize(WAHBitSet other) {
        int size = 0;

        // ensure that they have the same bit length.
        if (this.numBits() < other.numBits()) {
            this.setBit(other.numBits() - 1, 0);
        } else if (this.numBits() > other.numBits()) {
            other.setBit(numBits() - 1, 0);
        }

        if (vec.size() > 0) {
            run xrun = new run(vec), yrun = new run(other.vec);
            xrun.decode();
            yrun.decode();
            do {
                if (xrun.nWords == 0) {
                    xrun.inc();
                    xrun.decode();
                }

                if (yrun.nWords == 0) {
                    yrun.inc();
                    yrun.decode();
                }

                if (xrun.isFill()) {
                    if (yrun.isFill()) {
                        int nWords = Math.min(xrun.nWords, yrun.nWords);
                        if ((xrun.fillWord & yrun.fillWord) == 1) {
                            size += nWords * MAXBITS;
                        }
                        xrun.nWords -= nWords;
                        yrun.nWords -= nWords;
                    } else {
                        size += countInRun(xrun, yrun);
                    }
                } else if (yrun.isFill()) {
                    size += countInRun(yrun, xrun);
                } else {
                    int val = xrun.get() & yrun.get();
                    if (val > 0) {
                        size += Integer.bitCount(val);
                    }
                    yrun.nWords = 0;
                    xrun.nWords = 0;
                }

            } while (!(xrun.end() && yrun.end()));
        }
        size += Integer.bitCount(this.active.val & other.active.val);

        return size;
    }

    private int numBits() {
        return ((nbits != 0 ? nbits : (nbits = doCount())) + active.nbits);
    }

    private int doCount() {
        nset = 0;
        nbits = 0;

        for (int i = 0; i < vec.size(); i++) {
            int v = vec.getQuick(i);
            if (!isAFill(v)) {
                nbits += MAXBITS;
                nset += Integer.bitCount(v);
            } else {
                int tmp = (v & MAXCNT) * MAXBITS;
                nbits += tmp;
                nset += tmp * (isOneFill(v) ? 1 : 0);
            }
        }

        return nbits;
    }

    private void setBit(int ind, int val) {
        assert val == 0 || val == 1;

        if (ind >= numBits()) {
            int diff = ind - numBits() + 1;
            if (active.nbits > 0) {
                if (ind + 1 >= nbits + MAXBITS) {
                    diff -= MAXBITS - active.nbits;
                    active.val <<= (MAXBITS - active.nbits);
                    if (diff == 0) {
                        active.val += (val != 0 ? 1 : 0);
                    }
                    appendLiteral();
                } else {
                    active.nbits += diff;
                    active.val <<= diff;
                    active.val += (val != 0 ? 1 : 0);
                    diff = 0;
                }
            }
            if (diff != 0) {
                int w = diff / MAXBITS;
                diff -= w * MAXBITS;
                if (diff != 0) {
                    if (w > 1) {
                        appendCounter(0, w);
                    } else if (w != 0) {
                        appendLiteral();
                    }
                    active.nbits = diff;
                    active.val += (val != 0 ? 1 : 0);
                } else if (val != 0) {
                    if (w > 2) {
                        appendCounter(0, w - 1);
                    } else if (w == 2) {
                        appendLiteral();
                    }
                    active.val = 1;
                    appendLiteral();
                } else if (w > 1) {
                    appendCounter(0, w);
                } else if (w != 0) {
                    appendLiteral();
                }
            }

            if (numBits() != ind + 1) {
                logger.warn("Warning");
            }

            if (nset != 0) {
                nset += (val != 0 ? 1 : 0);
            }

            return;
        } else if (ind >= nbits) { // modify an active bit
            int u = active.val;
            if (val != 0) {
                active.val |= (1 << (active.nbits - (ind - nbits) - 1));
            } else {
                active.val &= ~(1 << (active.nbits - (ind - nbits) - 1));
            }
            if (nset != 0 && (u != active.val)) {
                nset += (val != 0 ? 1 : -1);
            }
            return;
        } else if (vec.size() * MAXBITS == nbits) { // uncompressed
            int i = ind / MAXBITS;
            int u = vec.get(i);
            int w = (1 << (SECONDBIT - (ind % MAXBITS)));

            if (val != 0) {
                vec.setQuick(i, u |= w);
            } else {
                vec.setQuick(i, u &= ~w);
            }
            if (nset != 0 && (vec.getQuick(i) != u)) {
                nset += (val != 0 ? 1 : -1);
            }
            return;
        }

        // the code after this has not been verified at all...
        // should proceed with caution.
        // compressed bit vector --
        // the bit to be modified is in vec
        if (RUN_UNTESTED_CODE) {
            int idx = 0;
            int compressed = 0, cnt = 0, ind1 = 0, ind0 = ind;
            int current = 0; // current bit value
            while ((ind0 > 0) && (idx < vec.size())) {
                int v = vec.getQuick(idx);

                if (isAFill(v)) { // a fill
                    cnt = ((v) & MAXCNT) * MAXBITS;
                    if (cnt > ind0) { // found the location
                        current = (isOneFill(v) ? 1 : 0);
                        compressed = 1;
                        ind1 = ind0;
                        ind0 = 0;
                    } else {
                        ind0 -= cnt;
                        ind1 = ind0;
                        ++idx;
                    }
                } else { // a literal word
                    cnt = MAXBITS;
                    if (MAXBITS > ind0) { // found the location
                        current = (1 & ((v) >>> (SECONDBIT - ind0)));
                        compressed = 0;
                        ind1 = ind0;
                        ind0 = 0;
                    } else {
                        ind0 -= MAXBITS;
                        ind1 = ind0;
                        ++idx;
                    }
                }
            } // while (ind...

            if (ind1 == 0) { // set current and compressed
                int v = vec.getQuick(idx);
                if (isAFill(v)) {
                    cnt = (v & MAXCNT) * MAXBITS;
                    current = (isOneFill(v) ? 1 : 0);
                    compressed = 1;
                } else {
                    cnt = MAXBITS;
                    current = (v >>> SECONDBIT);
                    compressed = 0;
                }
            }

            if (ind0 > 0) // has not found the right location yet.
            {
                if (ind0 < active.nbits) { // in the active word
                    ind1 = (1 << (active.nbits - ind0 - 1));
                    if (val != 0) {
                        active.val |= ind1;
                    } else {
                        active.val &= ~ind1;
                    }
                } else { // extends the current bit vector
                    ind1 = ind0 - active.nbits - 1;
                    appendWord(HEADER0 | (ind1 / MAXBITS));
                    for (ind1 %= MAXBITS; ind1 > 0; --ind1) {
                        addOneBit(0);
                    }
                    addOneBit(val != 0 ? 1 : 0);
                }
                if (nset != 0) {
                    nset += val != 0 ? 1 : -1;
                }
                return;
            }

            // locate the bit to be changed, lots of work hidden here
            if (current == val) {
                return; // nothing to do
            }
            int v = vec.getQuick(idx);
            // need to actually modify the bit
            if (compressed == 0) {
                // toggle a single bit of a literal word
                v ^= (1 << (SECONDBIT - ind1));
                vec.setQuick(idx, v);
            } else if (ind1 < MAXBITS) {
                // bit to be modified is in the first word, two pieces
                --v;
                vec.set(idx, v);
                if ((v & MAXCNT) == 1) {
                    v = (current != 0) ? ALLONES : 0;
                    vec.setQuick(idx, v);
                }
                int w = 1 << (SECONDBIT - ind1);
                if (val == 0) {
                    w ^= ALLONES;
                }

                vec.beforeInsert(idx, w);
                idx++;
            } else if (cnt - ind1 <= MAXBITS) {
                // bit to be modified is in the last word, two pieces
                --(v);
                vec.setQuick(idx, v);
                if ((v & MAXCNT) == 1) {
                    v = (current != 0) ? ALLONES : 0;
                    vec.setQuick(idx, v);
                }
                int w = 1 << (cnt - ind1 - 1);
                if (val == 0) {
                    w ^= ALLONES;
                }
                ++idx;
                vec.beforeInsert(idx, w);
            } else { // the counter breaks into three pieces
                int u[] = new int[2], w;
                u[0] = ind1 / MAXBITS;
                w = (v & MAXCNT) - u[0] - 1;
                u[1] = 1 << (SECONDBIT - ind1 + u[0] * MAXBITS);
                if (val == 0) {
                    u[0] = (u[0] > 1) ? (HEADER1 | u[0]) : (ALLONES);
                    u[1] ^= ALLONES;
                    w = (w > 1) ? (HEADER1 | w) : (ALLONES);
                } else {
                    u[0] = (u[0] > 1) ? (HEADER0 | u[0]) : 0;
                    w = (w > 1) ? (HEADER0 | w) : 0;
                }
                vec.setQuick(idx, w);
                vec.beforeInsertAllOf(idx, Arrays.asList(u));
            }

            if (nset != 0) {
                nset += val != 0 ? 1 : -1;
            }
        } else {
            throw new AssertionError("Untested code detected, would rather die than run this");
        }
    }

    private static void chewUpRun(run xrun, WAHBitSet ret, run yrun) {
        if (xrun.fillWord == 0) {
            // trye to advance, and see how much you can advance.
            int inc = yrun.inc(xrun.nWords);
            ret.appendFill(inc, 0);
            xrun.nWords -= inc;
        } else {
            while (true) {
                int v = yrun.get();
                if (isAFill(v)) {
                    yrun.decode();
                    return;
                }
                ret.appendCompressed(v);
                xrun.nWords--;
                if (xrun.nWords == 0) {
                    break;
                }
                yrun.inc();
            }
        }
        yrun.nWords = 0;
    }

    private static int countInRun(run xrun, run yrun) {
        int count = 0;
        if (xrun.fillWord == 0) {
            xrun.nWords -= yrun.inc(xrun.nWords);
        } else {
            while (true) {
                int v = yrun.get();
                if (isAFill(v)) {
                    yrun.decode();
                    return count;
                }
                if (v > 0) {
                    count += Integer.bitCount(v);
                }
                xrun.nWords--;
                if (xrun.nWords == 0) {
                    break;
                }
                yrun.inc();
            }
        }
        yrun.nWords = 0;
        return count;
    }

    private void appendCompressed(int v) {
        vec.add(v);
        nbits += MAXBITS;
        nset = 0;
    }

    private WAHBitSet genericOp(OpType op, WAHBitSet other) {
        WAHBitSet ret = new WAHBitSet();

        // ensure that they have the same bit length.
        if (this.numBits() < other.numBits()) {
            this.setBit(other.numBits() - 1, 0);
        } else if (this.numBits() > other.numBits()) {
            other.setBit(numBits() - 1, 0);
        }

        if (vec.size() > 0) {
            run xrun = new run(vec), yrun = new run(other.vec);
            xrun.decode();
            yrun.decode();
            do {
                if (xrun.nWords == 0) {
                    xrun.inc();
                    xrun.decode();
                }

                if (yrun.nWords == 0) {
                    yrun.inc();
                    yrun.decode();
                }

                if (xrun.isFill()) {
                    if (yrun.isFill()) {
                        int nWords = Math.min(xrun.nWords, yrun.nWords);
                        ret.appendFill(nWords, getOpResult(op, xrun.fillWord, yrun.fillWord));
                        xrun.nWords -= nWords;
                        yrun.nWords -= nWords;
                    } else {
                        ret.active.val = getOpResult(op, xrun.fillWord, yrun.get());
                        ret.appendLiteral();
                        --xrun.nWords;
                        yrun.nWords = 0;
                    }
                } else if (yrun.isFill()) {
                    ret.active.val = getOpResult(op, yrun.fillWord, xrun.get());
                    ret.appendLiteral();
                    yrun.nWords--;
                    xrun.nWords = 0;
                } else {
                    ret.active.val = getOpResult(op, xrun.get(), yrun.get());
                    ret.appendLiteral();
                    yrun.nWords = 0;
                    xrun.nWords = 0;
                }

            } while (!(xrun.end() && yrun.end()));
        }
        ret.active.val = getOpResult(op, this.active.val, other.active.val);
        ret.active.nbits = this.active.nbits;

        ret.doCount();
        return ret;
    }

    private void appendFill(int nWords, int v) {
        if (vec.isEmpty()) {
            if (v == 0) {
                vec.add(HEADER0 | nWords);
            } else {
                vec.add(HEADER1 | nWords);
            }
        } else if (nWords > 1) {
            int back = getBack();
            if (v == 0) {
                if (isZeroFill(back)) {
                    setBack(back + nWords);
                } else {
                    vec.add(HEADER0 | nWords);
                }
            } else if (isOneFill(back)) {
                setBack(back + nWords);
            } else {
                vec.add(HEADER1 | nWords);
            }
        } else {
            active.val = v != 0 ? ALLONES : 0;
            appendLiteral();
        }
    }

    private static int getOpResult(OpType op, int x, int y) {
        switch (op) {
            case And:
                return x & y;
            case AndNot:
                return x & ~y;
            case Or:
                return x | y;
        }

        throw new UnsupportedOperationException(op.toString());
    }

    private void addOneBit(int w) {
        if (RUN_UNTESTED_CODE) {
            int nb1, nb2;
            int cps = (w >>> MAXBITS);
            nset = 0;
            if (active.nbits != 0) { // active contains some uncompressed bits
                int w1;
                nb1 = active.nbits;
                nb2 = MAXBITS - active.nbits;
                active.val <<= nb2;
                if (cps != 0) { // incoming bits are comporessed
                    int b2 = (isOneFill(w) ? 1 : 0);
                    if (b2 != 0) {
                        w1 = (1 << nb2) - 1;
                        active.val |= w1;
                    }
                    appendLiteral();
                    nb2 = (w & MAXCNT) - 1;
                    if (nb2 > 1) {        // append a counter
                        appendCounter(b2, nb2);
                    } else if (nb2 == 1) {
                        if (b2 != 0) {
                            active.val = ALLONES;
                        }
                        appendLiteral();
                    }
                    active.nbits = nb1;
                    active.val = ((1 << nb1) - 1) * b2;
                } else { // incoming bits are not compressed
                    w1 = (w >>> nb1);
                    active.val |= w1;
                    appendLiteral();
                    w1 = (1 << nb1) - 1;
                    active.val = (w & w1);
                    active.nbits = nb1;
                }
            } // end of the case where there are active bits
            else if (cps != 0) { // no active bit
                int b2 = (isOneFill(w) ? 1 : 0);
                nb2 = (w & MAXCNT);
                if (nb2 > 1) {
                    appendCounter(b2, nb2);
                } else if (nb2 == 1) {
                    if (b2 != 0) {
                        active.val = ALLONES;
                    }
                    appendLiteral();
                }
            } else { // no active bits
                // new word is a raw bit pattern, simply add the word
                active.val = w;
                appendLiteral();
            }
        } else {
            throw new AssertionError("Untested code detected, would rather die than run this");
        }
    }

    private void appendWord(int w) {
        if (RUN_UNTESTED_CODE) {
            int nb1, nb2;
            int cps = (w >>> MAXBITS);
            nset = 0;
            if (active.nbits != 0) { // active contains some uncompressed bits
                int w1;
                nb1 = active.nbits;
                nb2 = MAXBITS - active.nbits;
                active.val <<= nb2;
                if (cps != 0) { // incoming bits are comporessed
                    boolean b2 = isOneFill(w);
                    if (b2) {
                        w1 = (1 << nb2) - 1;
                        active.val |= w1;
                    }
                    appendLiteral();
                    nb2 = (w & MAXCNT) - 1;
                    if (nb2 > 1) {        // append a counter
                        appendCounter(b2 ? 1 : 0, nb2);
                    } else if (nb2 == 1) {
                        if (b2) {
                            active.val = ALLONES;
                        }
                        appendLiteral();
                    }
                    active.nbits = nb1;
                    active.val = ((1 << nb1) - 1) * (b2 ? 1 : 0);
                } else { // incoming bits are not compressed
                    w1 = (w >>> nb1);
                    active.val |= w1;
                    appendLiteral();
                    w1 = (1 << nb1) - 1;
                    active.val = (w & w1);
                    active.nbits = nb1;
                }
            } // end of the case where there are active bits
            else if (cps != 0) { // no active bit
                int b2 = (isOneFill(w) ? 1 : 0);
                nb2 = (w & MAXCNT);
                if (nb2 > 1) {
                    appendCounter(b2, nb2);
                } else if (nb2 == 1) {
                    if (b2 != 0) {
                        active.val = ALLONES;
                    }
                    appendLiteral();
                }
            } else { // no active bits
                // new word is a raw bit pattern, simply add the word
                active.val = w;
                appendLiteral();
            }
        } else {
            throw new AssertionError("Untested code detected, would rather die than run this");
        }
    }

    private void appendLiteral() {
        if (vec.size() == 0) {
            vec.add(active.val);
        } else {
            int back = vec.getQuick(vec.size() - 1);

            if (active.val == 0) {// incoming word is zero
                if (back == 0) {
                    setBack(HEADER0 + 2);
                } else if (isZeroFill(back)) {
                    setBack(++back);
                } else {
                    vec.add(active.val);
                }
            } else if (active.val == ALLONES) {// incoming word is allones
                if (back == ALLONES) {
                    setBack((HEADER1 | 2));
                } else if (isOneFill(back)) {
                    setBack(++back);
                } else {
                    vec.add(active.val);
                }
            } else { // incoming word contains a mixture of bits
                vec.add(active.val);
            }
        }
        nbits += MAXBITS;
        active.reset();
        nset = 0;
    }

    private void appendCounter(int val, int cnt) {
        int head = 2 + val;
        int w = (head << SECONDBIT) + cnt;

        nbits += cnt * MAXBITS;
        if (vec.isEmpty()) {
            vec.add(w);
        } else {
            int back = getBack();
            if ((back >>> SECONDBIT) == head) {
                back += cnt;
                setBack(back);
            } else if ((back == ALLONES) && head == 3) {
                setBack(w + 1);
            } else if ((back == 0) && head == 2) {
                setBack(w + 1);
            } else {
                vec.add(w);
            }
        }
    }

    private int getBack() {
        return vec.getQuick(vec.size() - 1);
    }

    private void setBack(int val) {
        vec.setQuick(vec.size() - 1, val);
    }

    private static boolean isZeroFill(int v) {
        return (HEADER1 & v) == HEADER0;
    }

    private static boolean isAFill(int v) {
        return (v & HEADER0) == HEADER0;
    }

    private static boolean isOneFill(int v) {
        return (v & HEADER1) == HEADER1;
    }

    private static boolean moreThanInUnsigned(int val1, int val2) {
        if (val1 >= 0 && val2 >= 0) {
            return val1 > val2;
        } else if (val1 < 0 && val2 >= 0) {
            return true;
        } else if (val1 < 0 && val2 < 0) {
            return (val1 ^ HEADER0) > (val2 ^ HEADER0);
        } else {
            return false;
        }
    }

    private static int fillOnes(int v, int[] ret) {
        int count = 0;

        while (v != 0) {
            int i = Integer.numberOfLeadingZeros(v);
            ret[count++] = i - 1;
            v ^= 1 << (31 - i);
            // System.out.println("Integer.toHexString(i) = " + Integer.toHexString(v));
        }

        return count;
    }

    enum OpType {
        And,
        Or,
        AndNot,
    }

    // lets try our best for inlining
    // this class
    private static final class ActiveWord implements Serializable {

        int val;    // the value
        int nbits;    // total number of bits

        void reset() {
            val = 0;
            nbits = 0;
        }

        boolean is_full() {
            return (nbits >= MAXBITS);
        }

        void append(int b) { // append a single bit
            val <<= 1;
            nbits++;
            val += b;
        }
    }

    // this class is used for operations on the compressed
    private static final class run {

        int idx;
        int fillWord;
        int nWords;
        boolean fill;
        private IntArrayList vec;

        public run(IntArrayList l) {
            this.vec = l;
        }

        public void decode() {
            int v = vec.get(idx);
            if (WAHBitSet.isAFill(v)) {
                fillWord = (isOneFill(v) ? ALLONES : 0);
                nWords = v & MAXCNT;
                fill = true;
            } else {
                nWords = 1;
                fill = false;
            }
        }

        public boolean end() {
            return idx >= vec.size() - 1;
        }

        public void begin() {
            idx = 0;
        }

        public void inc() {
            idx++;
        }

        public int get() {
            return vec.getQuick(idx);
        }

        public boolean isFill() {
            return fill;
        }

        // this function tries to advance the run by nWords
        // sometimes it is not possible, in that case it returns how much it advanced.
        public int inc(int nWords) {
            int orig = nWords;
            nWords--;
            while (nWords > 0) {
                ++idx;
                int v = vec.get(idx);
                if (isAFill(v)) {
                    int words = v & MAXCNT;
                    if (words > nWords) {
                        --idx;
                        break;
                    }
                    nWords -= words;
                } else {
                    nWords--;
                }
            }

            return orig - nWords;
        }
    }

    /**
     * The IndexSet stores positions of bits that are one.
     *
     * It decodes one word of the bitvector at a time. For a fill of ones, the
     * function <tt>isRange</tt> returns true, otherwise it returns false. If
     * isRange returns true, the position of the first bit is pointed by the
     * pointer returned by function <tt>indices</tt>, and there are
     * <tt>nIndices</tt>
     * consecutive ones. If <tt>isRange</tt> returns false, there are
     * <tt>nIndices</tt>
     * bits that are one and the positions of these bits are stored in the array
     * returned by function <tt>indices</tt>.
     *
     * @see IndexSet#isRange
     * @see IndexSet#indices
     * @see IndexSet#nIndices
     */
    public class IndexSet {

        private int idx; // the crrent index in the vector.
        private int nind; // number of indices
        private int ind[] = new int[32]; // the returned array.
        boolean activeFilled = false; // just to check the last one

        public IndexSet() {
        }

        /**
         * @return true if the current word is a range word.
         */
        public boolean isRange() {
            return (nind >= MAXBITS);
        }

        /**
         * @return the indexes in the current word.
         */
        public int[] indices() {
            return ind;
        }

        /**
         * @return the number of indexes in the index array that are filled.
         */
        public int nIndices() {
            return nind;
        }

        /**
         * @return the raw current word.
         */
        private int currentWord() {
            return vec.get(idx);
        }

        public boolean hasMore() {
            return !activeFilled || idx <= vec.size() - 1;
        }

        /**
         * this function fills the indices and the nIndices.
         */
        public void next() {
            if (!hasMore()) { // already at the end
                nind = 0;
                return;
            }

            // the index of the next position
            int index0 = idx != 0 ? ((ind[0] + (nind > MAXBITS ? nind : MAXBITS)) / MAXBITS) * MAXBITS : 0;

            // reset the index count
            nind = 0;
            while (idx < vec.size()) {
                // get the current word.
                int v = currentWord();

                if (isOneFill(v)) { // 1-fill
                    // get the number of indexes.
                    nind = (v & MAXCNT) * MAXBITS;
                    ind[1] = index0 + nind;
                    ind[0] = index0;
                    ++idx;
                    return;
                } else if (isZeroFill(v)) { // 0-fill
                    // don't return, just skip to the next word.
                    index0 += (v & MAXCNT) * MAXBITS;
                    idx++;
                } else if (v > 0) { // non-zero literal word
                    if (v < ALLONES) { // a mixture of 0 and 1
                        nind = fillOnes(v, ind);

                        // add the base index.
                        for (int i = 0; i < nind; i++) {
                            ind[i] += index0;
                        }
                    } else { // all 1s
                        // ok, we have 31 1 bits.
                        nind = MAXBITS;
                        ind[0] = index0;
                        ind[1] = index0 + nind;
                    }

                    // update the next idx.
                    ++idx;
                    return;
                } else { // zero word
                    index0 += MAXBITS;
                    ++idx;
                }
            } // while (it < end)

            // deal with the active word
            if (active.nbits > 0 && active.val > 0) {
                // a non-zero active word
                // similar to the v > 0 case above.
                int j = (active.val << (MAXBITS - active.nbits));
                nind = fillOnes(j, ind);
                for (int i = 0; i < nind; i++) {
                    ind[i] += index0;
                }
            }

            // if you get here then you are done.
            activeFilled = true;
            // set the index accordingly.
            idx = vec.size() + 1;
            return;
        }
    }

    private class WAHIterator implements Iterator {

        private IndexSet is;
        private int max;
        private int current;
        private long[] pos;
        private long range_start;
        private long range_end;
        private int pos_idx;
        private boolean range;

        public WAHIterator() {
            is = new IndexSet();
            max = cardinality();
            current = 0;

            fill();
        }

        private void fill() {
            if (current < max) {
                is.next();

                int next[] = is.indices();
                range = is.isRange();
                if (range) {
                    range_start = next[0];
                    range_end = next[1];
                } else {
                    pos_idx = 0;
                    pos = new long[is.nIndices()];
                    for (int i = 0; i < pos.length; i++) {
                        pos[i] = next[i];
                    }
                }
            }
        }

        private int nextDoc() {
            current++;

            int ret = 0;
            if (range) {
                ret = (int) range_start++;
                if (range_start == range_end) {
                    fill();
                }
            } else {
                ret = (int) pos[pos_idx++];
                if (pos_idx == pos.length) {
                    fill();
                }
            }

            return ret;
        }

        public boolean hasNext() {
            return current < max;
        }

        public Object next() {
            return nextDoc();
        }

        public void remove() {
        }
    }

    public static void main(String[] args) {
        try {
            WAHBitSet set = new WAHBitSet();
            set.set(10);
            set.set(20);
            set.set(2000);
            set.set(20000);
            System.out.println(set.cardinality());
 
          //  assertEquals(set.cardinality(), 4);

            WAHBitSet set1 = new WAHBitSet();
            set1.set(1000);
            set1.set(20000);

           // assertTrue(set1.andSize(set) > 0);
            WAHBitSet and = set.and(set1);
            WAHBitSet or = set.or(set1);

           // assertEquals(and.cardinality(), 1);
           // assertEquals(or.cardinality(), 5);

            or.set(20003);
           // assertTrue(or.get(20000));
           // assertTrue(or.get(20003));

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

}
