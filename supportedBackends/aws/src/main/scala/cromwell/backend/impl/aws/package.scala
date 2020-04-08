package cromwell.backend.impl

import java.io.ByteArrayOutputStream
import java.util.zip.GZIPOutputStream

import cats.data.ReaderT
import com.google.common.io.BaseEncoding
import cromwell.cloudsupport.aws.auth.AwsAuthMode
import software.amazon.awssdk.awscore.client.builder.AwsClientBuilder
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.batch.model.KeyValuePair

import scala.language.higherKinds

package object aws {

  type Aws[F[_], A] = ReaderT[F, AwsBatchAttributes, A]

  def sanitize(name: String): String ={
    // Up to 128 letters (uppercase and lowercase), numbers, hyphens, and underscores are allowed.
    // We'll replace all invalid characters with an underscore
    name.replaceAll("[^A-Za-z0-9_\\-]", "_").slice(0, 128)
  }

  def buildKVPair(key: String, value: String): KeyValuePair =
    KeyValuePair.builder.name(key).value(value).build

  def gzipKeyValuePair(prefix: String, data: String): KeyValuePair = {
    // This limit is somewhat arbitrary. The idea here is that if there are
    // a reasonable amount of inputs, we'll just stick with "AWS_CROMWELL_INPUTS".
    // This allows for easier user debugging as everything is in plain text.
    // However, if there are a ton of inputs, then a) the user doesn't want
    // to wade through all that muck when debugging, and b) we want to provide
    // as much room for the rest of the container overrides (of which there
    // are currently none, but this function doesn't know that).
    //
    // This whole thing is setup because there is an undisclosed AWS Batch
    // limit of 8k for container override JSON when it hits the endpoint.
    // In the few circumstances where there are a ton of inputs, we can
    // hit this limit. This particular value is highly compressible, though,
    // so we get a ton of runway this way. Incidentally, the overhead of
    // gzip+base64 can be larger on very small inputs than the input data
    // itself, so we don't want it super-small either. The proxy container
    // is designed to handle either of these variables.
    val lim = 512 // This is completely arbitrary and I think 512 may even be
    // too big when looking at the value manually in the console
    data.length() match {
      case len if len <= lim => buildKVPair(prefix, data)
      case _ => buildKVPair(prefix + "_GZ", gzip(data))
    }
  }

  def gzip(data: String): String = {
    val byteArrayOutputStream = new ByteArrayOutputStream()
    val gzipOutputStream = new GZIPOutputStream(byteArrayOutputStream)
    gzipOutputStream.write(data.getBytes("UTF-8"))
    gzipOutputStream.close()

    BaseEncoding.base64().encode(byteArrayOutputStream.toByteArray())
  }

  /**
    * Generic method that, given a client builder, will configure that client with auth and region and return the client
    * to you
    * @param builder a builder for the client, the type of the builder will dictate the type of client
    * @param awsAuthMode an optional authorization mode
    * @param configRegion an optional region
    * @tparam BuilderT the type of builder (which dictates the type of client returned
    * @tparam ClientT the type of the client that you will get back
    * @return a configured client for the AWS service
    */
  def configureClient[BuilderT <: AwsClientBuilder[BuilderT, ClientT], ClientT](builder: AwsClientBuilder[BuilderT, ClientT],
                                                                                 awsAuthMode: Option[AwsAuthMode],
                                                                                 configRegion: Option[Region]): ClientT = {
    awsAuthMode.foreach { awsAuthMode =>
      builder.credentialsProvider(awsAuthMode.provider())
    }
    configRegion.foreach(builder.region)
    builder.build
  }
}
