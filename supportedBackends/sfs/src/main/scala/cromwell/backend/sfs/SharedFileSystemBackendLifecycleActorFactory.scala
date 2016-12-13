package cromwell.backend.sfs

import akka.actor.{ActorRef, Props}
import cats.data.Validated.{Invalid, Valid}
import cromwell.backend._
import cromwell.backend.standard.StandardLifecycleActorFactory
import cromwell.core.Dispatcher
import cromwell.core.Dispatcher._
import cromwell.core.path.{DefaultPathBuilderFactory, PathBuilderFactory}
import cromwell.filesystems.gcs.{GcsPathBuilderFactory, GoogleConfiguration}
import lenthall.exception.MessageAggregation
import net.ceedubs.ficus.Ficus._
import wdl4s.TaskCall
import wdl4s.expression.WdlStandardLibraryFunctions

/**
  * A factory that can be extended for any shared file system implementation.
  *
  * See the SharedFileSystemAsyncJobExecutionActor for more info.
  */
trait SharedFileSystemBackendLifecycleActorFactory extends StandardLifecycleActorFactory {

  override def jobIdKey: String = SharedFileSystemJob.JobIdKey

  /**
    * If the backend sets a gcs authentication mode, try to create a PathBuilderFactory with it.
    */
  lazy val gcsPathBuilderFactory: Option[GcsPathBuilderFactory] = {
    configurationDescriptor.backendConfig.as[Option[String]]("filesystems.gcs.auth") map { configAuth =>
      GoogleConfiguration(configurationDescriptor.globalConfig).auth(configAuth) match {
        case Valid(auth) => GcsPathBuilderFactory(auth)
        case Invalid(error) => throw new MessageAggregation {
          override def exceptionContext: String = "Failed to parse gcs auth configuration"

          override def errorMessages: Traversable[String] = error.toList
        }
      }
    }
  }

  lazy val pathBuilderFactories: List[PathBuilderFactory] = List(gcsPathBuilderFactory, Option(DefaultPathBuilderFactory)).flatten

  /**
    * Returns the initialization class, or by default uses the `SharedFileSystemInitializationActor`.
    *
    * @return the initialization class.
    */
  def initializationActorClass: Class[_ <: SharedFileSystemInitializationActor] =
    classOf[SharedFileSystemInitializationActor]

  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor, calls: Set[TaskCall],
                                                serviceRegistryActor: ActorRef): Option[Props] = {
    val params = SharedFileSystemInitializationActorParams(serviceRegistryActor, workflowDescriptor,
      configurationDescriptor, calls, pathBuilderFactories)
    Option(Props(initializationActorClass, params).withDispatcher(Dispatcher.BackendDispatcher))
  }

  override def cacheHitCopyingActorProps = Option(cacheHitCopyingActorInner _)

  def cacheHitCopyingActorInner(jobDescriptor: BackendJobDescriptor,
                                initializationDataOption: Option[BackendInitializationData],
                                serviceRegistryActor: ActorRef): Props = {
    Props(
      new SharedFileSystemCacheHitCopyingActor(
        jobDescriptor, configurationDescriptor, initializationDataOption, serviceRegistryActor)
    ).withDispatcher(BackendDispatcher)
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData]):
  WdlStandardLibraryFunctions = {
    SharedFileSystemExpressionFunctions(workflowDescriptor, configurationDescriptor, jobKey, initializationData)
  }
}
