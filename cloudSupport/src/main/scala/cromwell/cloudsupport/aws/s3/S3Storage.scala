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

import software.amazon.awssdk.services.s3.S3AdvancedConfiguration
import software.amazon.awssdk.services.s3.S3Client
import software.amazon.awssdk.core.auth.{AwsCredentials, StaticCredentialsProvider}
import com.typesafe.config.ConfigFactory
import net.ceedubs.ficus.Ficus._ // scalastyle:ignore

object S3Storage {
  val DefaultConfiguration = {
    val accelerateModeEnabled = ConfigFactory.load().as[Option[Boolean]]("s3.accelerate-mode").getOrElse(false)
    val dualstackEnabled = ConfigFactory.load().as[Option[Boolean]]("s3.dual-stack").getOrElse(false)
    val pathStyleAccessEnabled = ConfigFactory.load().as[Option[Boolean]]("s3.path-style-access").getOrElse(false)

    S3AdvancedConfiguration.builder
      .accelerateModeEnabled(accelerateModeEnabled)
      .dualstackEnabled(dualstackEnabled)
      .pathStyleAccessEnabled(pathStyleAccessEnabled)
      .build
  }

  def s3Client(configuration: S3AdvancedConfiguration, credentials: AwsCredentials): S3Client = {
    S3Client.builder
      .advancedConfiguration(configuration)
      .credentialsProvider(StaticCredentialsProvider.create(credentials))
      .build
  }

  def s3Client(credentials: AwsCredentials): S3Client = {
    s3Client(s3AdvancedConfiguration(), credentials)
  }

  def s3AdvancedConfiguration(accelerateModeEnabled: Boolean = false,
                              dualstackEnabled: Boolean = false,
                              pathStyleAccessEnabled: Boolean = false): S3AdvancedConfiguration = {

    S3AdvancedConfiguration.builder
      .accelerateModeEnabled(accelerateModeEnabled)
      .dualstackEnabled(dualstackEnabled)
      .pathStyleAccessEnabled(pathStyleAccessEnabled)
      .build
  }
}
