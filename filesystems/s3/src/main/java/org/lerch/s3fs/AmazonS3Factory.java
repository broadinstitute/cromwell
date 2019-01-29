package org.lerch.s3fs;

import software.amazon.awssdk.auth.credentials.AwsBasicCredentials;
import software.amazon.awssdk.auth.credentials.AwsCredentials;
import software.amazon.awssdk.auth.credentials.AwsCredentialsProvider;
import software.amazon.awssdk.auth.credentials.DefaultCredentialsProvider;
import software.amazon.awssdk.auth.credentials.StaticCredentialsProvider;
import software.amazon.awssdk.core.client.config.ClientOverrideConfiguration;
import software.amazon.awssdk.http.SdkHttpClient;
import software.amazon.awssdk.http.apache.ApacheHttpClient;
import software.amazon.awssdk.services.s3.S3Client;
import software.amazon.awssdk.services.s3.S3ClientBuilder;
import software.amazon.awssdk.services.s3.S3Configuration;

import java.net.URI;
import java.util.Properties;


/**
 * Factory base class to create a new AmazonS3 instance.
 */
public abstract class AmazonS3Factory {

    public static final String ACCESS_KEY = "s3fs_access_key";
    public static final String SECRET_KEY = "s3fs_secret_key";
    public static final String REQUEST_METRIC_COLLECTOR_CLASS = "s3fs_request_metric_collector_class";
    public static final String CONNECTION_TIMEOUT = "s3fs_connection_timeout";
    public static final String MAX_CONNECTIONS = "s3fs_max_connections";
    public static final String MAX_ERROR_RETRY = "s3fs_max_retry_error";
    public static final String PROTOCOL = "s3fs_protocol";
    public static final String PROXY_DOMAIN = "s3fs_proxy_domain";
    public static final String PROXY_HOST = "s3fs_proxy_host";
    public static final String PROXY_PASSWORD = "s3fs_proxy_password";
    public static final String PROXY_PORT = "s3fs_proxy_port";
    public static final String PROXY_USERNAME = "s3fs_proxy_username";
    public static final String PROXY_WORKSTATION = "s3fs_proxy_workstation";
    public static final String SOCKET_SEND_BUFFER_SIZE_HINT = "s3fs_socket_send_buffer_size_hint";
    public static final String SOCKET_RECEIVE_BUFFER_SIZE_HINT = "s3fs_socket_receive_buffer_size_hint";
    public static final String SOCKET_TIMEOUT = "s3fs_socket_timeout";
    public static final String USER_AGENT = "s3fs_user_agent";
    public static final String SIGNER_OVERRIDE = "s3fs_signer_override";
    public static final String PATH_STYLE_ACCESS = "s3fs_path_style_access";

    /**
     * Build a new Amazon S3 instance with the URI and the properties provided
     * @param uri URI mandatory
     * @param props Properties with the credentials and others options
     * @return S3Client
     */
    public S3Client getS3Client(URI uri, Properties props) {
        S3ClientBuilder builder = S3Client.builder();
        if (uri != null && uri.getHost() != null)
            builder.endpointOverride(uri);

        builder.credentialsProvider(getCredentialsProvider(props))
               .httpClient(getSdkHttpClient(props))
               .serviceConfiguration(getConfiguration(props))
               .overrideConfiguration(getOverrideConfiguration(props));
               //.region(getRegion(props));

        return createS3Client(builder);
    }

    /**
     * should return a new S3Client
     *
     * @return {@link software.amazon.awssdk.services.s3.S3Client}
     */
    protected abstract S3Client createS3Client(S3ClientBuilder builder);

    protected AwsCredentialsProvider getCredentialsProvider(Properties props) {
        AwsCredentialsProvider credentialsProvider;
        if (props.getProperty(ACCESS_KEY) == null && props.getProperty(SECRET_KEY) == null)
            credentialsProvider = DefaultCredentialsProvider.create();
        else
            credentialsProvider = StaticCredentialsProvider.create(getAWSCredentials(props));
        return credentialsProvider;
    }

    protected AwsCredentials getAWSCredentials(Properties props) {
        return AwsBasicCredentials.create(props.getProperty(ACCESS_KEY), props.getProperty(SECRET_KEY));
    }

    protected SdkHttpClient getSdkHttpClient(Properties props) {
        // TODO: custom http configuration based on properties
        return ApacheHttpClient.builder().build();
    }

    protected S3Configuration getConfiguration(Properties props) {
        // TODO: custom configuration based on properties
        return S3Configuration.builder().build();
    }

    protected ClientOverrideConfiguration getOverrideConfiguration(Properties props) {
        // TODO: custom configuration based on properties
        return ClientOverrideConfiguration.builder().build();
    }
}
