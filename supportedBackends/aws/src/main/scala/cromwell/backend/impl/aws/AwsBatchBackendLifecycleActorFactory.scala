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

import akka.actor.{ActorRef, Props}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor, JobExecutionMap}
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardFinalizationActor, StandardFinalizationActorParams, StandardInitializationActor, StandardInitializationActorParams, StandardLifecycleActorFactory}
import cromwell.core.CallOutputs
import wom.graph.CommandCallNode
import org.slf4j.LoggerFactory

case class AwsBatchBackendLifecycleActorFactory(
  name: String,
  configurationDescriptor: BackendConfigurationDescriptor)
    extends StandardLifecycleActorFactory {
  lazy val Log = LoggerFactory.getLogger(AwsBatchBackendLifecycleActorFactory.getClass)
  override lazy val initializationActorClass: Class[_ <: StandardInitializationActor]
    = classOf[AwsBatchInitializationActor]

  override lazy val asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor]
    = classOf[AwsBatchAsyncBackendJobExecutionActor]

  override lazy val finalizationActorClassOption: Option[Class[_ <: StandardFinalizationActor]]
    = Option(classOf[AwsBatchFinalizationActor])

  override lazy val jobIdKey: String
    = AwsBatchAsyncBackendJobExecutionActor.AwsBatchOperationIdKey

  val configuration = new AwsBatchConfiguration(configurationDescriptor)

  override def workflowInitializationActorParams(
                                                  workflowDescriptor: BackendWorkflowDescriptor,
                                                  ioActor: ActorRef,
                                                  calls: Set[CommandCallNode],
                                                  serviceRegistryActor: ActorRef,
                                                  restart: Boolean): StandardInitializationActorParams = {
    Log.debug("Initializing AwsBatchBackendLifecycleActorFactory")
    AwsBatchInitializationActorParams(workflowDescriptor, ioActor, calls, configuration, serviceRegistryActor, restart)
  }

  override def workflowFinalizationActorParams(
                                                workflowDescriptor: BackendWorkflowDescriptor,
                                                ioActor: ActorRef,
                                                calls: Set[CommandCallNode],
                                                jobExecutionMap: JobExecutionMap,
                                                workflowOutputs: CallOutputs,
                                                initializationDataOption: Option[BackendInitializationData]): StandardFinalizationActorParams = {
    // The `AwsBatchInitializationActor` will only return a non-`Empty`
    // `AwsBatchBackendInitializationData` from a successful `beforeAll`
    // invocation.  HOWEVER, the finalization actor is created regardless
    // of whether workflow initialization was successful or not.  So the
    // finalization actor must be able to handle an empty
    // `AwsBatchBackendInitializationData` option, and there is no `.get`
    // on the initialization data as there is with the execution or cache
    // hit copying actor methods.
    AwsBatchFinalizationActorParams(workflowDescriptor, ioActor, calls, configuration, jobExecutionMap, workflowOutputs, initializationDataOption)
  }

  override def backendSingletonActorProps(serviceRegistryActor: ActorRef): Option[Props] = {
    Some(AwsBatchSingletonActor.props(None))
  }
}
