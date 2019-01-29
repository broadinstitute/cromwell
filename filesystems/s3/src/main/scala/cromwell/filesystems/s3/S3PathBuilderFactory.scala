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
package cromwell.filesystems.s3

import akka.actor.ActorSystem
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cromwell.cloudsupport.aws.AwsConfiguration
import cromwell.cloudsupport.aws.auth.AwsAuthMode
import cromwell.cloudsupport.aws.s3.S3Storage
import cromwell.core.path.PathBuilderFactory
import cromwell.core.WorkflowOptions
import net.ceedubs.ficus.Ficus._

import scala.concurrent.{ExecutionContext, Future}
import software.amazon.awssdk.auth.credentials.AwsCredentials

// The constructor of this class is required to be Config, Config by cromwell
// So, we need to take this config and get the AuthMode out of it
final case class S3PathBuilderFactory private(globalConfig: Config, instanceConfig: Config)
  extends PathBuilderFactory {

  // Grab the authMode out of configuration
  val conf: AwsConfiguration = AwsConfiguration(globalConfig)
  val authModeAsString: String = instanceConfig.as[String]("auth")
  val authModeValidation: ErrorOr[AwsAuthMode] = conf.auth(authModeAsString)
  val authMode = authModeValidation.unsafe(s"Failed to get authentication mode for $authModeAsString")

  def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[S3PathBuilder] = {
    S3PathBuilder.fromAuthMode(authMode, S3Storage.DefaultConfiguration,  options, conf.region)
  }

  // Ignores the authMode and creates an S3PathBuilder using the passed credentials directly.
  // Can be used when the Credentials are already available.
  def fromCredentials(options: WorkflowOptions, credentials: AwsCredentials): S3PathBuilder = {
    S3PathBuilder.fromCredentials(credentials, S3Storage.DefaultConfiguration, options, conf.region)
  }
}

object S3PathBuilderFactory {
  def apply(globalConfig: Config, instanceConfig: Config): S3PathBuilderFactory = {
    new S3PathBuilderFactory(globalConfig, instanceConfig)
  }
}
