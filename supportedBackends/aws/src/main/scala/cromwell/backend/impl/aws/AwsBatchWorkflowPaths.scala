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

package cromwell.backend.impl.aws

import software.amazon.awssdk.auth.credentials.AwsCredentials
import akka.actor.ActorSystem
import com.typesafe.config.Config
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.cloudsupport.aws.s3.S3Storage
import cromwell.core.WorkflowOptions
import cromwell.core.path.{Path, PathBuilder}
import cromwell.filesystems.s3.S3PathBuilder

import scala.language.postfixOps

object AwsBatchWorkflowPaths {
  private val RootOptionKey = "aws_s3_root"
  private val AuthFilePathOptionKey = "auth_bucket"
}

case class AwsBatchWorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor,
                            credentials: AwsCredentials,
                            configuration: AwsBatchConfiguration)(implicit actorSystem: ActorSystem) extends WorkflowPaths {

  override lazy val executionRootString: String =
    workflowDescriptor.workflowOptions.getOrElse(AwsBatchWorkflowPaths.RootOptionKey, configuration.root)

  private val workflowOptions: WorkflowOptions = workflowDescriptor.workflowOptions

  private val pathBuilder: S3PathBuilder = configuration.pathBuilderFactory.fromCredentials(workflowOptions, credentials)

  val authFilePath: Path = {
    // The default auth file bucket is always at the root of the root workflow
    val defaultBucket = executionRoot.resolve(workflowDescriptor.rootWorkflow.name).resolve(workflowDescriptor.rootWorkflowId.toString)
    val bucket = workflowDescriptor.workflowOptions.get(AwsBatchWorkflowPaths.AuthFilePathOptionKey) getOrElse defaultBucket.pathAsString

    val pathBuilderWithAuth = S3PathBuilder.fromCredentials(
      credentials,
      S3Storage.DefaultConfiguration,
      workflowOptions,
      configuration.awsConfig.region
    )

    val authBucket = pathBuilderWithAuth.build(bucket) recover {
      case ex => throw new Exception(s"Invalid s3 auth_bucket path $bucket", ex)
    } get

    authBucket.resolve(s"${workflowDescriptor.rootWorkflowId}_auth.json")
  }

  override def toJobPaths(workflowPaths: WorkflowPaths, jobKey: BackendJobDescriptorKey): AwsBatchJobPaths = {
    new AwsBatchJobPaths(workflowPaths.asInstanceOf[AwsBatchWorkflowPaths], jobKey)
  }

  override protected def withDescriptor(workflowDescriptor: BackendWorkflowDescriptor): WorkflowPaths = this.copy(workflowDescriptor = workflowDescriptor)

  override def config: Config = configuration.configurationDescriptor.backendConfig
  override def pathBuilders: List[PathBuilder] = List(pathBuilder)
}
