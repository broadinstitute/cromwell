package cromwell.backend.google.pipelines.batch

import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardLifecycleActorFactory}
//import akka.actor.{ActorRef, Props}
//import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor, JobExecutionMap}
//import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardFinalizationActor, StandardFinalizationActorParams, StandardInitializationActor, StandardInitializationActorParams, StandardLifecycleActorFactory}
//import cromwell.core.CallOutputs


/**
  * https://github.com/broadinstitute/cromwell/blob/develop/supportedBackends/aws/src/main/scala/cromwell/backend/impl/aws/AwsBatchAsyncBackendJobExecutionActor.scala
  * https://github.com/Burwood/cromwell/blob/dspeck-11-14-2022/backend/src/main/scala/cromwell/backend/async/AsyncBackendJobExecutionActor.scala
  */


/**
  * Factory to create `Actor` objects to manage the lifecycle of a backend job on AWS Batch. This factory provides an
  * object from the `AwsBatchAsyncBackendJobExecutionActor` class to create and manage the job.
  * @param name Factory name
  * @param configurationDescriptor configuration descriptor for the backend
  */

case class GcpBatchBackendLifecycleFactory(name: String,  configurationDescriptor: BackendConfigurationDescriptor) extends StandardLifecycleActorFactory {

  val gcpBatchConfig = new GcpBatchConfiguration(configurationDescriptor)
  println(gcpBatchConfig)

  /**
    * Returns the key to use for storing and looking up the job id.
    *
    * @return the key to use for storing and looking up the job id.
    */
  override def jobIdKey: String = "gcp_batch"

  /**
    * Returns the asynchronous executor class.
    *
    * @return the asynchronous executor class.
    */
  override def asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] = classOf[GcpBatchAsyncBackendJobExecutionActor]
}
