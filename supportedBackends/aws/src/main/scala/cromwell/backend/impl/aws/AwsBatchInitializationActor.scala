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

import java.io.IOException

import akka.actor.ActorRef
import software.amazon.awssdk.auth.credentials.AwsCredentials
import cromwell.filesystems.s3.batch.S3BatchCommandBuilder
import cromwell.backend.standard.{StandardInitializationActor,
                                  StandardInitializationActorParams,
                                  StandardValidatedRuntimeAttributesBuilder}
import cromwell.backend.{BackendConfigurationDescriptor,
                         BackendWorkflowDescriptor}
import cromwell.core.io.DefaultIoCommandBuilder
import cromwell.core.io.AsyncIoActorClient
import cromwell.core.path.Path
import wom.graph.CommandCallNode
import scala.concurrent.Future

case class AwsBatchInitializationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  ioActor: ActorRef,
  calls: Set[CommandCallNode],
  configuration: AwsBatchConfiguration,
  serviceRegistryActor: ActorRef,
  restarting: Boolean
) extends StandardInitializationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = configuration.configurationDescriptor
}

object AwsBatchInitializationActor {
  private case class AuthFileAlreadyExistsException(path: Path) extends IOException(s"Failed to upload authentication file at $path:" +
    s" there was already a file at the same location and this workflow was not being restarted.")
}

class AwsBatchInitializationActor(params: AwsBatchInitializationActorParams)
  extends StandardInitializationActor(params) with AsyncIoActorClient {

  override lazy val ioActor = params.ioActor
  private val configuration = params.configuration
  override implicit val system = context.system

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    AwsBatchRuntimeAttributes.runtimeAttributesBuilder(configuration)

  private lazy val credentials: Future[AwsCredentials] =
    Future { configuration.awsAuth.credential(_ => "") }

  override lazy val workflowPaths: Future[AwsBatchWorkflowPaths] = for {
    creds <- credentials
  } yield new AwsBatchWorkflowPaths(workflowDescriptor, creds, configuration)

  override lazy val initializationData: Future[AwsBatchBackendInitializationData] = for {
    workflowPaths <- workflowPaths
    creds <- credentials
  } yield AwsBatchBackendInitializationData(workflowPaths, runtimeAttributesBuilder, configuration, creds)

  override lazy val ioCommandBuilder =  {
    val conf = Option(configuration) match {
      case Some(cf) => cf
      case None =>  new  AwsBatchConfiguration(params.configurationDescriptor)
    }
    conf.fileSystem match {
      case  AWSBatchStorageSystems.s3 =>  S3BatchCommandBuilder
      case _ =>   DefaultIoCommandBuilder
    }
  }
}
