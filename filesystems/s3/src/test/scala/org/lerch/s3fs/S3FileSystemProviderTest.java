package org.lerch.s3fs;


import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.mockito.InOrder;
import org.mockito.Mockito;
import software.amazon.awssdk.services.s3.S3Client;
import software.amazon.awssdk.services.s3.model.*;

import java.io.IOException;
import java.net.URI;
import java.nio.file.StandardCopyOption;
import java.time.Instant;
import java.util.Properties;

import static org.junit.Assert.*;
import static org.mockito.Mockito.*;

/**
 * An integration test that confirms that the S3FileSystemProvider makes appropriate calls to an s3Client which
 * for test purposes is mocked. Also tests that responses from the client are appropriately handled. The pattern can
 * be repeated and extended to increase test coverage of the s3fs package.
 */
public class S3FileSystemProviderTest {

    S3FileSystemProvider s3fsProvider;
    S3FileSystem s3fs;
    S3Client s3Client;


    @Before
    public void setUp() throws Exception {
        s3Client = Mockito.mock(S3Client.class);


        s3fsProvider = new S3FileSystemProvider();
        s3fs = s3fsProvider.createFileSystem(new URI("s3.amazonaws.com"), new Properties(), s3Client);
    }

    @After
    public void tearDown() throws Exception {
        if(s3fs != null) {
            s3fsProvider.close(s3fs);
        }
    }

    /**
     * Configures the standard responses that the mock client should make when requests are made
     */
    private void standardMockSetup(){
        when(s3Client.headObject(any(HeadObjectRequest.class))).thenReturn(
                HeadObjectResponse.builder()
                        .lastModified(Instant.now())
                        .eTag("fake-etag")
                        .contentLength(10L)
                        .storageClass("standard")
                        .build());

        when(s3Client.getObjectAcl(any(GetObjectAclRequest.class))).thenReturn(
                GetObjectAclResponse.builder()
                        .owner(builder -> builder.displayName("alice").id("id"))
                        .build());
    }

    /**
     * Setup the mock s3 client to contain  a 6GB object and respond to multipart upload requests
     */
    private void largeObjectMockSetup(){
        Long sixGB = 1024L * 1024L * 1024L * 6L;
        when(s3Client.headObject(any(HeadObjectRequest.class))).thenReturn(
                HeadObjectResponse.builder()
                        .lastModified(Instant.now())
                        .eTag("fake-etag")
                        .contentLength(sixGB)
                        .storageClass("standard")
                        .build());

        when(s3Client.getObjectAcl(any(GetObjectAclRequest.class))).thenReturn(
                GetObjectAclResponse.builder()
                        .owner(builder -> builder.displayName("alice").id("id"))
                        .build());

        when(s3Client.createMultipartUpload(any(CreateMultipartUploadRequest.class))).thenReturn(
                CreateMultipartUploadResponse.builder()
                        .uploadId("id")
                        .build());

        when(s3Client.uploadPartCopy(any(UploadPartCopyRequest.class))).thenReturn(
                UploadPartCopyResponse.builder()
                        .copyPartResult(builder -> builder.eTag("fake-etag").lastModified(Instant.now()))
                        .build());
    }

    @Test
    public void getScheme() {
        assertEquals("Scheme should be 's3'", "s3", s3fsProvider.getScheme());
    }

    @Test
    public void delete() throws IOException {
        standardMockSetup();

        S3Path path = s3fs.getPath("/testbucket", "/file/name");
        s3fsProvider.delete(path);

        //ensure that delete call to the fs results in a delete call to the client
        verify(s3Client, atLeastOnce()).deleteObject(any(DeleteObjectRequest.class));
    }

    @Test
    public void copyLargeObject() throws IOException {
        largeObjectMockSetup();
        s3fsProvider.copy(
                s3fs.getPath("/testbucket", "/file/name"),
                s3fs.getPath("/testbucket", "/file/name2"),
                StandardCopyOption.REPLACE_EXISTING);

        InOrder inOrder = Mockito.inOrder(s3Client);
        inOrder.verify(s3Client, atLeastOnce()).createMultipartUpload(any(CreateMultipartUploadRequest.class));
        inOrder.verify(s3Client, atLeastOnce()).uploadPartCopy(any(UploadPartCopyRequest.class));
        inOrder.verify(s3Client, atLeastOnce()).completeMultipartUpload(any(CompleteMultipartUploadRequest.class));

        //no more than 10K parts are allowed
        verify(s3Client, atMost(10000)).uploadPartCopy(any(UploadPartCopyRequest.class));
    }

    @Test
    public void copy() throws IOException{
        standardMockSetup();
        s3fsProvider.copy(
                s3fs.getPath("/testbucket", "/file/name"),
                s3fs.getPath("/testbucket", "/file/name2"),
                StandardCopyOption.REPLACE_EXISTING);

        verify(s3Client, atLeastOnce()).copyObject(any(CopyObjectRequest.class));

    }

}