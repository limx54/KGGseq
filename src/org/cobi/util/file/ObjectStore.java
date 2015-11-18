/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.file;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.LinkedList;
import java.util.List;

/**
 *
 * @author modified from
 */
public class ObjectStore {

    private File file;
    private FileInputStream fis;
    private FileOutputStream fos;
    private FileChannel channel;
    private final static int _DELETE = -1;
    private final static int BSIZE = 1024;
    private final static int BSIZE_FLAG = 4;
    private long NEW_POSITION;

    public ObjectStore(String fileName) throws IOException {
        this.file = new File(fileName);
    }

    /**
     * Creates a channel to read from a file
     */
    public void read() throws IOException {
        FileInputStream fis = new FileInputStream(file);
        this.channel = fis.getChannel();
    }

    /**
     * Creates a channel to write to a file in appended mode
     */
    public void write() throws IOException {
        FileOutputStream fos = new FileOutputStream(file, true);
        this.channel = fos.getChannel();
    }

    /**
     * Returns the object at the given index.
     */
    public Object getObject(int index) throws IOException,
            ClassNotFoundException {
        seekToIndex(index);
        return readObject();
    }

    public void close() throws IOException {
        this.channel.close();
    }

    /**
     * Delete the object at index from the store.
     */
    public void deleteObject(int index) throws IOException {
        seekToIndex(index);

        NEW_POSITION = (long) (this.channel.position() + BSIZE_FLAG);
        this.channel.position(NEW_POSITION);

        writeDeleteFlag(_DELETE);
    }

    /**
     * Adds this object to the store. Always added at the end.
     */
    public void addObject(Object o) throws IOException {
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        ObjectOutputStream oos = new ObjectOutputStream(bos);
        oos.writeObject(o);
        byte[] buffer = bos.toByteArray();

// move to the last position
        this.channel.position(this.channel.size());

        writeSizeInfo(buffer);

        writeDeleteFlag(0);

// Write the actual object as byteBuffer
        this.channel.write(ByteBuffer.wrap(buffer));

    }

    private void writeDeleteFlag(int flag) throws IOException {
// write an int flag to identify, whether deleted or not, if delete, its
// _DELETE, i.e -1 else 0.
        ByteBuffer flagBuffer = ByteBuffer.allocateDirect(BSIZE_FLAG);
        flagBuffer.asIntBuffer().put(flag);
        this.channel.write(flagBuffer);
    }

    private void writeSizeInfo(byte[] buffer) throws IOException {
// write the number of bytes in buffer to be written
        ByteBuffer sizeBuffer = ByteBuffer.allocate(BSIZE_FLAG);
        sizeBuffer.asIntBuffer().put(buffer.length);
        this.channel.write(sizeBuffer);
    }

    /**
     * Returns the number of objects currently in the store.
     */
    public int getStoredCount() throws IOException {

        this.channel.position(0);
        int counter = 0;
        while (this.channel.position() < this.channel.size()) {

            int skipBytes = readSizeInfo();
            int deleteFlag = readDeleteFlag();

            if (deleteFlag != _DELETE) {
                counter++;
            }

// skip number of bytes
            NEW_POSITION = (long) (this.channel.position() + skipBytes);
            this.channel.position(NEW_POSITION);
        }
        System.out.println("Total Objects Written to a file is: " + counter);
        return counter;
    }

    /**
     * @return @throws IOException
     */
    private int readDeleteFlag() throws IOException {
// read flag information
        ByteBuffer flagBuffer = ByteBuffer.allocate(BSIZE_FLAG);
        this.channel.read(flagBuffer);
        flagBuffer.rewind();
        int deleteFlag = flagBuffer.getInt();
        return deleteFlag;
    }

    /**
     * @return @throws IOException
     */
    private int readSizeInfo() throws IOException {
// read and skip size information
        ByteBuffer sizeBuffer = ByteBuffer.allocate(BSIZE_FLAG);
        this.channel.read(sizeBuffer);
        sizeBuffer.rewind();
        int skipBytes = sizeBuffer.getInt();
        return skipBytes;
    }

    /**
     * Returns all the objects in the store as a List.
     */
    public List getObjects() throws IOException, ClassNotFoundException {
        List list = new LinkedList();
        this.channel.position(0);
        while (this.channel.position() < this.channel.size()) {

            int skipBytes = readSizeInfo();
            int deleteFlag = readDeleteFlag();

            if (deleteFlag != _DELETE) {
                this.channel.position(this.channel.position() - 8);
                list.add(readObject());
            }

// skip number of bytes
            NEW_POSITION = (long) (this.channel.position() + skipBytes);
            this.channel.position(NEW_POSITION);

        }
        return list;
    }

    /**
     * Return al the objects inside specific range mentioned in the store as a
     * List
     */
    public List getObjects(int fromIndex, int toIndex) throws IOException,
            ClassNotFoundException {
        List list = new LinkedList();
        for (int i = fromIndex + 1; i < toIndex; i++) {
            list.add(getObject(i));
        }
        return list;
    }

    /**
     * Compacts the underlying file by removing all the dead space aka deleted
     * records.
     */
    public void compact() throws IOException {
// TO DO
    }

    private void seekToIndex(int index) throws IOException {
        this.channel.position(0);
        int counter = 0;
        /*
         * seek to the correct index skipping those we don't want. deleted
         * records do not count as space
         */
        while (counter < index) {
            if (this.channel.position() == this.channel.size()) {
                throw new IOException("Object " + index + " not present.");
            }

            int skipBytes = readSizeInfo();
            int deleteFlag = readDeleteFlag();

// skip number of bytes
            NEW_POSITION = (long) (this.channel.position() + skipBytes);
            this.channel.position(NEW_POSITION);

            if (deleteFlag != _DELETE) {
                counter++;
            }

        }
        /*
         * now in theory we are in the right place but this index could be
         * deleted. keep reading until we are at an existing record
         */
        boolean found = false;
        while (!found) {

            if (this.channel.position() == this.channel.size()) {
                throw new IOException("Object " + index + " not present.");
            }

            int skipBytes = readSizeInfo();
            int deleteFlag = readDeleteFlag();

            if (deleteFlag != _DELETE) {
                this.channel.position(this.channel.position() - 8);
                found = true;
            } else {
// skip number of bytes
                NEW_POSITION = (long) (this.channel.position() + skipBytes);
                this.channel.position(NEW_POSITION);
            }

        }
    }

// does the actual byte to Object comversion for getObject(int) and
// getObjects()
    private Object readObject() throws IOException, ClassNotFoundException {
        int bytes = readSizeInfo();
        readDeleteFlag();
        byte[] buff = new byte[bytes];
        this.channel.read(ByteBuffer.wrap(buff));
        ObjectInputStream ois = new ObjectInputStream(new ByteArrayInputStream(buff));
        return ois.readObject();
    }
}
