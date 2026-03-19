/*
 * Copyright 2018 Amazon.com, Inc. or its affiliates.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 *  BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 *  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 *  THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 *  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 *  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 *  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */
package cromwell.cloudsupport.aws.s3

import com.typesafe.config.ConfigFactory
import java.net.URI
import net.ceedubs.ficus.Ficus._
import scala.annotation.nowarn
import software.amazon.awssdk.auth.credentials.AwsCredentialsProvider
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.s3.{S3Client, S3Configuration}

object S3Storage {
  val DefaultConfiguration = {
    val accelerateModeEnabled = ConfigFactory.load().as[Option[Boolean]]("s3.accelerate-mode").getOrElse(false)
    val dualstackEnabled = ConfigFactory.load().as[Option[Boolean]]("s3.dual-stack").getOrElse(false)
    val pathStyleAccessEnabled = ConfigFactory.load().as[Option[Boolean]]("s3.path-style-access").getOrElse(false)

    @nowarn("msg=method dualstackEnabled in trait Builder is deprecated")
    val builder = S3Configuration.builder
      .accelerateModeEnabled(accelerateModeEnabled)
      .dualstackEnabled(dualstackEnabled)
      .pathStyleAccessEnabled(pathStyleAccessEnabled)

    builder.build
  }

  /**
   * Build an S3Client.
   *
   * @param endpointUri optional custom S3-compatible endpoint (e.g. OVH, MinIO).
   *                    When set, path-style access is forced automatically because:
   *                    (a) non-AWS services do not support virtual-hosted-style bucket addressing,
   *                    (b) AWS SDK v2 would otherwise prepend the bucket name to the custom hostname.
   */
  def s3Client(
    configuration: S3Configuration,
    provider: AwsCredentialsProvider,
    region: Option[Region],
    endpointUri: Option[URI] = None
  ): S3Client = {
    val builder = S3Client.builder.credentialsProvider(provider)
    // For a custom endpoint, rebuild the S3Configuration with path-style access enabled.
    // Non-AWS S3-compatible services require path-style (https://endpoint/bucket/key)
    // rather than virtual-hosted-style (https://bucket.endpoint/key).
    val finalConfig = endpointUri match {
      case Some(_) => configuration.toBuilder.pathStyleAccessEnabled(true).build()
      case None    => configuration
    }
    builder.serviceConfiguration(finalConfig)
    region.foreach(builder.region)
    endpointUri.foreach(builder.endpointOverride)
    builder.build
  }

  def s3Client(provider: AwsCredentialsProvider, region: Option[Region]): S3Client =
    s3Client(s3Configuration(), provider, region)

  def s3Configuration(accelerateModeEnabled: Boolean = false,
                      dualstackEnabled: Boolean = false,
                      pathStyleAccessEnabled: Boolean = false
  ): S3Configuration = {

    @nowarn("msg=method dualstackEnabled in trait Builder is deprecated")
    val builder = S3Configuration.builder
      .accelerateModeEnabled(accelerateModeEnabled)
      .dualstackEnabled(dualstackEnabled)
      .pathStyleAccessEnabled(pathStyleAccessEnabled)

    builder.build()
  }
}
