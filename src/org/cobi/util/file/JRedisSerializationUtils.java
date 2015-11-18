/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.file;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import org.nustaq.serialization.FSTConfiguration;

/**
 *
 * @author mxli
 */
public class JRedisSerializationUtils {

    public JRedisSerializationUtils() {
    }

    static FSTConfiguration configuration = FSTConfiguration // .createDefaultConfiguration();  
            .createStructConfiguration();

    public static byte[] serialize(Object obj) {
        return configuration.asByteArray(obj);
    }

    public static Object unserialize(byte[] sec) {
        return configuration.asObject(sec);
    }

    public static byte[] kryoSerizlize(Object obj) {
        Kryo kryo = new Kryo();
        byte[] buffer = new byte[2048];
        try (
                Output output = new Output(buffer);)
        {

            kryo.writeClassAndObject(output, obj);
            return output.toBytes();
        } catch (Exception e) {
        }
        return buffer;
    }

    static Kryo kryo = new Kryo();

    public static Object kryoUnSerizlize(byte[] src) {
        try (
                Input input = new Input(src);) {
            return kryo.readClassAndObject(input);
        } catch (Exception e) {
        }
        return kryo;
    }

    // jdk原生序列换方案  
    public static byte[] jdkserialize(Object obj) {
        try (ByteArrayOutputStream baos = new ByteArrayOutputStream();
                ObjectOutputStream oos = new ObjectOutputStream(baos);) {
            oos.writeObject(obj);
            return baos.toByteArray();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static Object jdkdeserialize(byte[] bits) {
        try (ByteArrayInputStream bais = new ByteArrayInputStream(bits);
                ObjectInputStream ois = new ObjectInputStream(bais);) {
            return ois.readObject();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
}
