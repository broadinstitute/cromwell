package org.lerch.s3fs.util;

import junit.framework.TestCase;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.mockito.InOrder;
import org.mockito.Mock;
import static org.mockito.Mockito.*;

import org.mockito.Spy;
import org.mockito.junit.MockitoJUnitRunner;
import software.amazon.awssdk.awscore.exception.AwsErrorDetails;
import software.amazon.awssdk.http.SdkHttpResponse;
import software.amazon.awssdk.services.s3.S3Client;
import software.amazon.awssdk.services.s3.model.GetBucketLocationResponse;
import software.amazon.awssdk.services.s3.model.HeadBucketResponse;
import software.amazon.awssdk.services.s3.model.S3Exception;

import java.util.NoSuchElementException;
import java.util.function.Consumer;

@RunWith(MockitoJUnitRunner.class)
public class S3ClientStoreTest extends TestCase {

    S3ClientStore instance;

    @Mock
    S3Client mockClient;

    @Spy
    final S3ClientStore spyInstance = S3ClientStore.getInstance();

    @Before
    public void setUp() throws Exception {
        super.setUp();
        instance = S3ClientStore.getInstance();
    }

    @Test
    public void testGetInstanceReturnsSingleton() {
        assertSame(S3ClientStore.getInstance(), instance);
    }

    @Test
    public void testGetClientForNullBucketName() {
        assertEquals(S3ClientStore.DEFAULT_CLIENT, instance.getClientForBucketName(null));
    }

    @Test
    public void testGetClientForEmptyBucketName() {
        assertEquals(S3ClientStore.DEFAULT_CLIENT, instance.getClientForBucketName(""));
        assertEquals(S3ClientStore.DEFAULT_CLIENT, instance.getClientForBucketName(" "));
    }

    @Test
    public void testGenerateClientWithNoErrors() {
        when(mockClient.getBucketLocation(any(Consumer.class)))
                .thenReturn(GetBucketLocationResponse.builder().locationConstraint("us-west-2").build());
        final S3Client s3Client = instance.generateClient("test-bucket", mockClient);
        assertNotNull(s3Client);

        ;
    }

    @Test
    public void testGenerateClientWith403Response() {
        // when you get a forbidden response from getBucketLocation
        when(mockClient.getBucketLocation(any(Consumer.class))).thenThrow(
                S3Exception.builder().statusCode(403).build()
        );
        // you should fall back to a head bucket attempt
        when(mockClient.headBucket(any(Consumer.class)))
                .thenReturn((HeadBucketResponse) HeadBucketResponse.builder()
                        .sdkHttpResponse(SdkHttpResponse.builder()
                                .putHeader("x-amz-bucket-region", "us-west-2")
                                .build())
                        .build());

        // which should get you a client
        final S3Client s3Client = instance.generateClient("test-bucket", mockClient);
        assertNotNull(s3Client);

        final InOrder inOrder = inOrder(mockClient);
        inOrder.verify(mockClient).getBucketLocation(any(Consumer.class));
        inOrder.verify(mockClient).headBucket(any(Consumer.class));
        inOrder.verifyNoMoreInteractions();
    }

    @Test
    public void testGenerateClientWith403Then301Responses(){
        // when you get a forbidden response from getBucketLocation
        when(mockClient.getBucketLocation(any(Consumer.class))).thenThrow(
                S3Exception.builder().statusCode(403).build()
        );
        // and you get a 301 response on headBucket
        when(mockClient.headBucket(any(Consumer.class))).thenThrow(
                S3Exception.builder()
                        .statusCode(301)
                        .awsErrorDetails(AwsErrorDetails.builder()
                                .sdkHttpResponse(SdkHttpResponse.builder()
                                        .putHeader("x-amz-bucket-region", "us-west-2")
                                        .build())
                                .build())
                        .build()
        );

        // then you should be able to get a client as long as the error response header contains the region
        final S3Client s3Client = instance.generateClient("test-bucket", mockClient);
        assertNotNull(s3Client);

        final InOrder inOrder = inOrder(mockClient);
        inOrder.verify(mockClient).getBucketLocation(any(Consumer.class));
        inOrder.verify(mockClient).headBucket(any(Consumer.class));
        inOrder.verifyNoMoreInteractions();
    }

    @Test
    public void testGenerateClientWith403Then301ResponsesNoHeader(){
        // when you get a forbidden response from getBucketLocation
        when(mockClient.getBucketLocation(any(Consumer.class))).thenThrow(
                S3Exception.builder().statusCode(403).build()
        );
        // and you get a 301 response on headBucket but no header for region
        when(mockClient.headBucket(any(Consumer.class))).thenThrow(
                S3Exception.builder()
                        .statusCode(301)
                        .awsErrorDetails(AwsErrorDetails.builder()
                                .sdkHttpResponse(SdkHttpResponse.builder()
                                        .build())
                                .build())
                        .build()
        );

        // then you should get a NoSuchElement exception when you try to get the header
        try {
            instance.generateClient("test-bucket", mockClient);
        } catch (Exception e) {
            assertEquals(NoSuchElementException.class, e.getClass());
        }

        final InOrder inOrder = inOrder(mockClient);
        inOrder.verify(mockClient).getBucketLocation(any(Consumer.class));
        inOrder.verify(mockClient).headBucket(any(Consumer.class));
        inOrder.verifyNoMoreInteractions();
    }

    @Test
    public void testCaching() {
        S3Client client = S3Client.create();
        doReturn(client).when(spyInstance).generateClient("test-bucket");

        final S3Client client1 = spyInstance.getClientForBucketName("test-bucket");
        verify(spyInstance).generateClient("test-bucket");
        assertSame(client1, client);

        S3Client differentClient = S3Client.create();
        assertNotSame(client, differentClient);

        lenient().doReturn(differentClient).when(spyInstance).generateClient("test-bucket");
        final S3Client client2 = spyInstance.getClientForBucketName("test-bucket");
        // same instance because second is cached.
        assertSame(client1, client2);
        assertSame(client2, client);
        assertNotSame(client2, differentClient);
    }
}
