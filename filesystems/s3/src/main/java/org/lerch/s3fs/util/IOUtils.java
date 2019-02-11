package org.lerch.s3fs.util;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * Utilities for streams
 */
public abstract class IOUtils {
    /**
     * get the stream content and return as a byte array
     *
     * @param is InputStream
     * @return byte array
     * @throws IOException if the stream is closed
     */
    public static byte[] toByteArray(InputStream is) throws IOException {
        ByteArrayOutputStream buffer = new ByteArrayOutputStream();

        int nRead;
        byte[] data = new byte[16384];

        while ((nRead = is.read(data, 0, data.length)) != -1) {
            buffer.write(data, 0, nRead);
        }

        buffer.flush();

        return buffer.toByteArray();
    }
}
