/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.text;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;
import static org.cobi.util.text.BGZFInputStream.BUF_SIZE;

/**
 *
 * @author mxli
 */
public class BZPartReader   {

    int spiderID;
    int intFormat=0;
    long longRemainSize = -1;
    long longEnd = -1;
    InputStream inputStream;
    LimitInputStream cis;

    int intRead = -1;
    byte[] bytBuffer = new byte[BUF_SIZE];
    int intLineStart = 0;
    int intCurrPos = 0;
    int intTempBufferLengthACC = 0;
    byte[] bytTempBuffer = new byte[BUF_SIZE];
    byte[] bytStartLine = null;
    byte[] bytEndLine = null;
    byte bytDelimiter = 9;
    int[] intSep;
    byte[] buf = new byte[BUF_SIZE];

    File input;

    public BZPartReader(File input, int spiderID, int intFormat, long longStart, long longEnd, char chrDelimiter) throws IOException {
        this.input = input;
        longRemainSize = (longEnd - longStart);
        this.spiderID = spiderID;
        this.intFormat = intFormat;

        if (intFormat == 0) {
            InputStream is = new BufferedInputStream(new FileInputStream(input));
            cis = new LimitInputStream(is, longEnd);
            cis.skip(longStart);
            inputStream = new GZIPInputStream(cis);
        } else {
            InputStream is = new BufferedInputStream(new FileInputStream(input));
            cis = new LimitInputStream(is, longEnd);
            cis.skip(longStart);
            inputStream = cis;
        }

        if (this.spiderID != 0) {
            bytStartLine = this.readLine();
            // System.out.println(new String(bytStartLine));
        }
        bytDelimiter = (byte) chrDelimiter;
    }

      public BZPartReader(File input,int spiderID, int intFormat, long longStart, long longEnd, char chrDelimiter, boolean boolRFL) throws IOException {
        longRemainSize = (longEnd - longStart);
        this.spiderID = spiderID;
        this.intFormat = intFormat;

        if (intFormat == 0) {
            InputStream is = new BufferedInputStream(new FileInputStream(input), 100000);
            cis = new LimitInputStream(is, longEnd);
            cis.skip(longStart);
            inputStream = new GZIPInputStream(cis);
        } else {
            InputStream is = new BufferedInputStream(new FileInputStream(input));
            cis = new LimitInputStream(is, longEnd);
            cis.skip(longStart);
            inputStream = cis;
        }

        if (this.spiderID == 0 && boolRFL) {//boolRFL<=>longStart!=0
            bytStartLine = this.readLine();
        }

        if (this.spiderID != 0) {
            bytStartLine = this.readLine();         
        }
        bytDelimiter = (byte) chrDelimiter;
    }

    public BZPartReader(File input,int spiderID, long longStart, long longEnd, int[] intSep, boolean hasHead) throws IOException {
        this.longEnd = longEnd;
        longRemainSize = (longEnd - longStart);
        this.spiderID = spiderID;
        this.intSep = intSep;

        InputStream is = new BufferedInputStream(new FileInputStream(input));
        cis = new LimitInputStream(is, longEnd);
        cis.skip(longStart);
        cis.mark((int) longRemainSize);
        if (intFormat == 0) {
            inputStream = new GZIPInputStream(cis);
        } else {
            inputStream = cis;
        }

        if (this.spiderID != 0) {
            this.readLine(intSep);
            bytStartLine = readLine(intSep);

        } else {
            if (hasHead) {
                this.readLine(intSep);
            }
            bytStartLine = readLine(intSep);
        }
        cis.reset();
    }

    public void closeInputStream() throws IOException {
        inputStream.close();
    }

    public synchronized byte[] readLine(int[] intDelPos) throws IOException {
        int intDelPosMarker = 0;
        intDelPos[intDelPosMarker] = 0;
        int len = intDelPos.length;
        do {
            if (intRead == -1) {
                intRead = (inputStream.read(bytBuffer));
                if (intRead == -1) {
                    //The end of the block is not a complete line. 
                    if (intTempBufferLengthACC != 0) {
                        bytBuffer = new byte[intTempBufferLengthACC + (bytEndLine == null ? 0 : bytEndLine.length) + 1];
                        System.arraycopy(bytTempBuffer, 0, bytBuffer, 0, intTempBufferLengthACC);
                        if (bytEndLine != null) {
                            System.arraycopy(this.bytEndLine, 0, bytBuffer, intTempBufferLengthACC, this.bytEndLine.length);
                        }
                        bytBuffer[this.bytBuffer.length - 1] = (byte) '\n';
                        intTempBufferLengthACC = 0;
                        intDelPosMarker = 0;
                    } else {
                        return null;
                    }
                }
            }

            intLineStart = intCurrPos;
            while (intCurrPos != intRead) {
                if (bytBuffer[intCurrPos] == bytDelimiter) {
                    ++intDelPosMarker;
                    if (intDelPosMarker < len) {
                        intDelPos[intDelPosMarker] = intCurrPos - intLineStart + intTempBufferLengthACC;
                    }
                } else if (bytBuffer[intCurrPos] == 10) {
                    //parse the line. 
                    int intLineLength = intCurrPos - intLineStart;//don't contaion \n
                    byte[] bytLine = null;
                    if (intTempBufferLengthACC != 0) {
                        // bytLine = new byte[intTempBufferLengthACC + intLineLength];
                        //System.arraycopy(bytTempBuffer, 0, bytLine, 0, intTempBufferLengthACC);                        
                        bytLine = Arrays.copyOfRange(bytTempBuffer, 0, intTempBufferLengthACC + intLineLength);
                        System.arraycopy(bytBuffer, intLineStart, bytLine, intTempBufferLengthACC, intLineLength);
                        intTempBufferLengthACC = 0;
                    } else {
                        // bytLine = new byte[intLineLength];
                        // System.arraycopy(bytBuffer, intLineStart, bytLine, 0, intLineLength);
                        bytLine = Arrays.copyOfRange(bytBuffer, intLineStart, intCurrPos);
                    }
                    intCurrPos++;
                    if (intCurrPos == intRead) {
                        intRead = -1;
                        intCurrPos = 0;
                    }

                    ++intDelPosMarker;
                    if (intDelPosMarker < len) {
                        //return 13 new line 10
                        //the new line or return line symbol is not included
                        if (intCurrPos > 1 && bytBuffer[intCurrPos - 2] == 13) {
                            intDelPos[intDelPosMarker] = bytLine.length - 1;
                        } else {
                            intDelPos[intDelPosMarker] = bytLine.length;
                        }
                        intDelPos[0] = intDelPosMarker;
                    } else {
                        intDelPos[0] = len;
                    }

                    return bytLine;
                }
                intCurrPos++;
            }

            //The buffer ends with imcomplete line. 
            int intTempBufferLength = intCurrPos - intLineStart;
            if ((bytTempBuffer.length - intTempBufferLengthACC) < intTempBufferLength) {
                bytTempBuffer = Arrays.copyOf(bytTempBuffer, bytTempBuffer.length + intTempBufferLength * 2);
            }
            System.arraycopy(bytBuffer, intLineStart, bytTempBuffer, intTempBufferLengthACC, intTempBufferLength);
            intTempBufferLengthACC += intTempBufferLength;
            intRead = -1;
            intCurrPos = 0;
        } while (intTempBufferLengthACC != 0);

        return null;
    }

    public synchronized byte[] readLine() throws IOException {
        do {
            if (intRead == -1) {
                intRead = (inputStream.read(bytBuffer));
                if (intRead == -1) {
                    //The end of the block is not a complete line. 
                    if (intTempBufferLengthACC != 0) {
                        bytBuffer = new byte[intTempBufferLengthACC + (bytEndLine == null ? 0 : bytEndLine.length) + 1];
                        System.arraycopy(bytTempBuffer, 0, bytBuffer, 0, intTempBufferLengthACC);
                        if (bytEndLine != null) {
                            System.arraycopy(this.bytEndLine, 0, bytBuffer, intTempBufferLengthACC, this.bytEndLine.length);
                        }
                        bytBuffer[this.bytBuffer.length - 1] = (byte) '\n';
                        intTempBufferLengthACC = 0;
                    } else {
                        return null;
                    }
                }
            }

            intLineStart = intCurrPos;
            while (intCurrPos != intRead) {
                if (bytBuffer[intCurrPos] == 10) {
                    //parse the line. 
                    int intLineLength = intCurrPos - intLineStart;//don't contaion \n
                    byte[] bytLine = null;
                    if (intTempBufferLengthACC != 0) {
                        // bytLine = new byte[intTempBufferLengthACC + intLineLength];
                        //System.arraycopy(bytTempBuffer, 0, bytLine, 0, intTempBufferLengthACC);                        
                        bytLine = Arrays.copyOfRange(bytTempBuffer, 0, intTempBufferLengthACC + intLineLength);
                        System.arraycopy(bytBuffer, intLineStart, bytLine, intTempBufferLengthACC, intLineLength);
                        intTempBufferLengthACC = 0;
                    } else {
                        // bytLine = new byte[intLineLength];
                        // System.arraycopy(bytBuffer, intLineStart, bytLine, 0, intLineLength);
                        bytLine = Arrays.copyOfRange(bytBuffer, intLineStart, intCurrPos);
                    }
                    intCurrPos++;
                    if (intCurrPos == intRead) {
                        intRead = -1;
                        intCurrPos = 0;
                    }
                    return bytLine;
                }
                intCurrPos++;
            }

            //The buffer ends with imcomplete line. 
            int intTempBufferLength = intCurrPos - intLineStart;
            if ((bytTempBuffer.length - intTempBufferLengthACC) < intTempBufferLength) {
                bytTempBuffer = Arrays.copyOf(bytTempBuffer, bytTempBuffer.length + intTempBufferLength * 2);
            }
            System.arraycopy(bytBuffer, intLineStart, bytTempBuffer, intTempBufferLengthACC, intTempBufferLength);
            intTempBufferLengthACC += intTempBufferLength;
            intRead = -1;
            intCurrPos = 0;
        } while (intTempBufferLengthACC != 0);

        return null;
    }

    public byte[] getStartLine() {
        return bytStartLine;
    }

    public long getNextBlock() throws IOException {
        long longNext = 0;
        //cis.skip(longCurr + 1);
        boolean boolNoFound = true;
        int n_buf = -1;
        int intContent = 0;
        do {
            intContent = cis.read(buf);
            if (intContent == -1) {
                return -1;
            }
            n_buf++;
            for (int id = 0; id < intContent - 4; id++) {
                if (buf[id] == 31 && buf[id + 1] == -117 && buf[id + 2] == 8 && buf[id + 3] == 4 && buf[id + 4] == 0) { //This should be used unsigned number or others. 
                    longNext += (id + n_buf * buf.length + 1);
                    boolNoFound = false;
                    break;
                }
            }
        } while (boolNoFound);
        return longNext;
    }

    protected boolean jumpTo(long longDict, int[] intSepTmp) throws IOException {
        longRemainSize -= longDict;
        // System.out.println(longRemainSize);
        if (longRemainSize <= 0) {
            bytStartLine = null;
            return false;
        }
        while (longDict > 0) {
            longDict -= cis.skip(longDict);
        }
        cis.mark((int) longRemainSize);
//            cis.skip(longDict);

        if (intFormat == 0) {
//            inputStream.close();
//            InputStream inputStream2=new GZIPInputStream(cis);
//            inputStream.close();

            inputStream = new GZIPInputStream(cis);
//            inputStream=inputStream2;
//            inputStream2.close();
        } else {
            inputStream = cis;
        }

        intRead = -1;
        intTempBufferLengthACC = 0;
        intCurrPos = 0;
        readLine(intSepTmp);
        bytStartLine = readLine(intSepTmp);

        // System.out.println(new String(bytStartLine));
        cis.reset();
//            inputStream.reset();
//            System.out.println("jump right!");
        return true;
    }

    public boolean prepareRead(long longDict) throws IOException {
        longRemainSize -= longDict;
        // System.out.println(longRemainSize);
        if (longRemainSize <= 28) {
            bytStartLine = null;
            return false;
        }
        while (longDict > 0) {
            longDict -= cis.skip(longDict);
//                longDict-=inputStream.skip(longDict);
        }
        cis.mark((int) longRemainSize);
        cis.skip(longDict);
//            inputStream.mark((int) longRemainSize);

        if (intFormat == 0) {
            // inputStream = new GZIPInputStream(cis);
        } else {
            //inputStream = cis;
        }

        intRead = -1;
        intTempBufferLengthACC = 0;
        intCurrPos = 0;

        return true;
    }
}
