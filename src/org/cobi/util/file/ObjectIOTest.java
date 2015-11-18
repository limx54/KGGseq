/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.file;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.RandomAccessFile;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import org.cobi.kggseq.entity.Variant;
import org.objenesis.strategy.StdInstantiatorStrategy;

/**
 *
 * @author mxli
 */
public class ObjectIOTest {

    public Object testReadBuffered(String fileName) throws Exception {
        Object test = null;
        ObjectInputStream objectInputStream = null;
        FileInputStream fos = null;
        BufferedInputStream bos = null;
        try {
            fos = new FileInputStream(fileName);
            bos = new BufferedInputStream(fos);
            objectInputStream = new ObjectInputStream(bos);
            test = objectInputStream.readObject();
        } finally {
            if (objectInputStream != null) {
                objectInputStream.close();
            }
            bos.close();
            fos.close();
        }
        return test;
    }
//source http://java.dzone.com/articles/fast-java-file-serialization

    public void testWriteBuffered(Object test, String fileName) throws IOException {
        ObjectOutputStream objectOutputStream = null;
        try {
            FileOutputStream fos = new FileOutputStream(fileName);
            BufferedOutputStream bos = new BufferedOutputStream(fos);
            objectOutputStream = new ObjectOutputStream(bos);
            objectOutputStream.writeObject(test);
            bos.close();
            fos.close();
        } finally {
            if (objectOutputStream != null) {
                objectOutputStream.close();
            }
        }
    }

    public void testWriteBuffered1(Object test, String fileName) throws IOException {
        ObjectOutputStream objectOutputStream = null;
        try {
            RandomAccessFile raf = new RandomAccessFile(fileName, "rw");
            FileOutputStream fos = new FileOutputStream(raf.getFD());
            objectOutputStream = new ObjectOutputStream(fos);
            objectOutputStream.writeObject(test);
            fos.close();
            raf.close();
        } finally {
            if (objectOutputStream != null) {
                objectOutputStream.close();
            }
        }
    }
    //some good packages http://blog.csdn.net/blognkliming/article/details/16333865

    static String inFile = "kggseq1/Chromosome.2.obj.0";
    //  static String inFile = "test.txt.gz";
    static String outFile = inFile + "tmp";

    //Conculstion: I spent one day in comparing performance of io and nio. I found io is much faster than nio on Windows
    //Probaly io has been optimized a lot for text parsing.
    public static void main(String[] args) throws Exception {
        ObjectIOTest oit = new ObjectIOTest();
        // oit.sytaxTest();

        Kryo kryo = new Kryo();
        Long startTime = System.nanoTime();
        Input input = new Input(new FileInputStream(inFile), 1024 * 1024);
        List<Variant> varList = new ArrayList<Variant>();
        kryo.setReferences(false);
        kryo.setRegistrationRequired(false);
        kryo.setInstantiatorStrategy(new StdInstantiatorStrategy());
        //kryo.setInstantiatorStrategy(new SerializingInstantiatorStrategy());
        
        kryo.register(ArrayList.class);
        kryo.register(Variant.class);
        Variant var = null;
        while (!input.eof()) {
            var = (Variant) kryo.readObject(input, Variant.class);
            if (var == null) {
                break;
            }
            varList.add(var);
            // varList.add(var);
        }
        System.out.println(varList.size());
        input.close();

        Long endTime = System.nanoTime();
        System.out.println("用时：" + (endTime - startTime) / 1000000000.0);
        startTime = System.nanoTime();

        endTime = System.nanoTime();
        System.out.println("用时：" + (endTime - startTime) / 1000000000.0);

    }

    static class Variant1 implements Serializable {

        public int refStartPosition;
        private String label;
        private String refGeneAnnot;
        private String gEncodeAnnot;
        private String knownGeneAnnot;
        private String ensemblGeneAnnot;
        private String refAllele;
        private String[] altAlleles;
        public float[] scores;
        public String geneSymb;
        public boolean isIndel = false;
        //-1 denotes this SNP does not exist in db; NA means db has this variant but no frequency information
        public float altAF = -1;
        public List<String> featureValues;
        public boolean isIBS = false;
        public byte smallestFeatureID = 17;//by default
        public int genotypeIndex = -1;
        private int affectedRefHomGtyNum = 0;
        private int affectedHetGtyNum = 0;
        private int affectedAltHomGtyNum = 0;
        private int unaffectedRefHomGtyNum = 0;
        private int unaffectedHetGtyNum = 0;
        private int unaffectedAltHomGtyNum = 0;

        public int getAffectedAltHomGtyNum() {
            return affectedAltHomGtyNum;
        }

        public void setAffectedAltHomGtyNum(int affectedAltHomGtyNum) {
            this.affectedAltHomGtyNum = affectedAltHomGtyNum;
        }

        public int getAffectedHetGtyNum() {
            return affectedHetGtyNum;
        }

        public void setAffectedHetGtyNum(int affectedHetGtyNum) {
            this.affectedHetGtyNum = affectedHetGtyNum;
        }

        public int getAffectedRefHomGtyNum() {
            return affectedRefHomGtyNum;
        }

        public void setAffectedRefHomGtyNum(int affectedRefHomGtyNum) {
            this.affectedRefHomGtyNum = affectedRefHomGtyNum;
        }

        public int getUnaffectedAltHomGtyNum() {
            return unaffectedAltHomGtyNum;
        }

        public void setUnaffectedAltHomGtyNum(int unaffectedAltHomGtyNum) {
            this.unaffectedAltHomGtyNum = unaffectedAltHomGtyNum;
        }

        public int getUnaffectedHetGtyNum() {
            return unaffectedHetGtyNum;
        }

        public void setUnaffectedHetGtyNum(int unaffectedHetGtyNum) {
            this.unaffectedHetGtyNum = unaffectedHetGtyNum;
        }

        public int getUnaffectedRefHomGtyNum() {
            return unaffectedRefHomGtyNum;
        }

        public void setUnaffectedRefHomGtyNum(int unaffectedRefHomGtyNum) {
            this.unaffectedRefHomGtyNum = unaffectedRefHomGtyNum;
        }

        public String getEnsemblGeneAnnot() {
            return ensemblGeneAnnot;
        }

        public void setEnsemblGeneAnnot(String ensemblGeneAnnot) {
            this.ensemblGeneAnnot = ensemblGeneAnnot;
        }

        public List<String> getFeatureValues() {
            return featureValues;
        }

        public int getAlleleEndPostion(int alleleIndex) throws Exception {
            if (!isIndel) {
                return refStartPosition;
            } else {
                if (alleleIndex >= altAlleles.length) {
                    throw new Exception("Allele index is out of boundery " + altAlleles.length);
                }
                String alta = altAlleles[alleleIndex];
                if (alta.startsWith("+")) {
                    return refStartPosition + refAllele.length() - 1;
                } else {
                    //deletion
                    return refStartPosition + refAllele.length() - 1;
                }
            }
        }

        public String getKnownGeneAnnot() {
            return knownGeneAnnot;
        }

        public void setKnownGeneAnnot(String knownGeneAnnot) {
            this.knownGeneAnnot = knownGeneAnnot;
        }

        public String getgEncodeAnnot() {
            return gEncodeAnnot;
        }

        public void setgEncodeAnnot(String gEncodeAnnot) {
            this.gEncodeAnnot = gEncodeAnnot;
        }

        public String getRefGeneAnnot() {
            return refGeneAnnot;
        }

        public void setRefGeneAnnot(String refGeneAnnot) {
            this.refGeneAnnot = refGeneAnnot;
        }

        //for some indels the refStartPosition is acturally not actually mutant allele position when the reference allele have multiple bases
        //this function report one base before the insertion or deletion allele
        public int getAlleleStartPostion(int alleleIndex) {
            if (!isIndel) {
                return refStartPosition;
            } else {
                if (alleleIndex >= altAlleles.length) {
                    //  throw new Exception("Allele index is out of boundery " + altAlleles.length);
                    return -1;
                }
                String alta = altAlleles[alleleIndex];
                if (alta.startsWith("+")) {
                    return refStartPosition + refAllele.length() - 1;
                } else {
                    //deletion
                    int index = alta.indexOf('-');
                    return refStartPosition + index - 1;
                }
            }
        }

        public String getAltAllele(int alleleIndex) {
            String alta = altAlleles[alleleIndex];
            if (!isIndel) {
                return alta;
            } else {
                if (alleleIndex >= altAlleles.length) {
                    //  throw new Exception("Allele index is out of boundery " + altAlleles.length);
                    return null;
                }
                alta = altAlleles[alleleIndex];
                if (alta.startsWith("+")) {
                    return alta.substring(refAllele.length());
                } else if (alta.startsWith("-")) {
                    //deletion
                    int index = alta.indexOf('-');
                    return refAllele.substring(index);
                } else {
                    StringBuilder sb = new StringBuilder();
                    for (int i = 0; i < alta.length(); i++) {
                        if (refAllele.charAt(i) != alta.charAt(i)) {
                            sb.append(alta.charAt(i));
                        }
                    }
                    return sb.toString();
                }
            }
        }

        public void setIsIndel(boolean isIndel) {
            this.isIndel = isIndel;
        }

        public void addFeatureValue(String val) {
            featureValues.add(val);
        }

        public Variant1() {
        }

        public Variant1(int physicalPosition, String refAllele, String[] altAlleles) {
            this.refStartPosition = physicalPosition;
            this.refAllele = refAllele;
            this.altAlleles = altAlleles;
            if (altAlleles != null) {
                for (String s : altAlleles) {
                    if (s.length() > 1) {
                        isIndel = true;
                        break;
                    }
                }
            }
            featureValues = new ArrayList<String>();
        }

        public void setFeatureValues(List<String> featureValues) {
            this.featureValues = featureValues;
        }

        public String[] getAltAlleles() {
            return altAlleles;
        }

        public void setAltAlleles(String[] altAlleles) {
            this.altAlleles = altAlleles;
        }

        public String getLabel() {
            return label;
        }

        public void setLabel(String label) {
            this.label = label;
        }

        public String getRefAllele() {
            return refAllele;
        }

        public void setRefAllele(String refAllele) {
            this.refAllele = refAllele;
        }
    }

    public void sytaxTest() throws Exception {
        Variant1 var = new Variant1(100, "12", new String[]{"we"});
        Kryo kryo = new Kryo();
        // ByteBuffer buffer = ByteBuffer.allocateDirect(512);
        var.refStartPosition = 1000;

        List<String> myList = new ArrayList<String>();
        myList.add("apples");
        myList.add("bananas");
        var.setFeatureValues(myList);
        Output output = new Output();
        byte[] buffer = new byte[1024 * 1024];
        output.setBuffer(buffer);
        output.setOutputStream(new FileOutputStream("test.obj"));
        kryo.register(Variant1.class);
        kryo.register(ArrayList.class);    //for arrayList used inside the bean
        kryo.writeClassAndObject(output, var);
        kryo.writeClassAndObject(output, var);
        kryo.writeClassAndObject(output, var);
        output.close();
        RandomAccessFile raf = new RandomAccessFile("test.obj", "r");
        Input input = new Input(new FileInputStream(raf.getFD()), 1024 * 1024);

        //  kryo.setReferences(false);
        //  kryo.setRegistrationRequired(false);
        //     kryo.setInstantiatorStrategy(new StdInstantiatorStrategy());
        // kryo.register(List.class);
        while (!input.eof()) {
            var = (Variant1) kryo.readObject(input, Variant1.class);
            if (var == null) {
                break;
            }
            System.out.println(var.refStartPosition);
            // varList.add(var);

        }
        //Variant1 newBean = kryo.readObject(input, Variant1.class);
        System.out.println(var.refStartPosition);
        System.out.println(var.refStartPosition);
    }
}
