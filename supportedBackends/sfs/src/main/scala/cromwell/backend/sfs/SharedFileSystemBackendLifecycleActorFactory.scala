package cromwell.backend.sfs

import cromwell.backend._
import cromwell.backend.standard._
import wdl4s.expression.WdlStandardLibraryFunctions

/**
  * A factory that can be extended for any shared file system implementation.
  *
  * See the SharedFileSystemAsyncJobExecutionActor for more info.
  */
trait SharedFileSystemBackendLifecycleActorFactory extends StandardLifecycleActorFactory {

  override def jobIdKey: String = SharedFileSystemAsyncJobExecutionActor.JobIdKey

  override lazy val standardCacheHitCopyingActorOption: Option[Class[_ <: StandardCacheHitCopyingActor]] = {
    Option(classOf[SharedFileSystemCacheHitCopyingActor])
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData]):
  WdlStandardLibraryFunctions = {
    SharedFileSystemExpressionFunctions(workflowDescriptor, configurationDescriptor, jobKey, initializationData)
  }
}
