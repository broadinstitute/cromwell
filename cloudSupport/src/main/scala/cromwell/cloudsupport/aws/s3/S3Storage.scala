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

  def s3Client(configuration: S3Configuration, provider: AwsCredentialsProvider, region: Option[Region]): S3Client = {
    val builder = S3Client.builder
      .serviceConfiguration(configuration)
      .credentialsProvider(provider)
    region.foreach(builder.region)
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
