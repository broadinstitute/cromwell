package cromwell.backend.standard

import akka.actor.{ActorRef, Props}
import cromwell.backend._
import cromwell.core.Dispatcher

/**
  * May be extended for using the standard sync/async backend pattern.
  */
trait StandardLifecycleActorFactory extends BackendLifecycleActorFactory {
  /**
    * Config values for the backend, and a pointer to the global config.
    *
    * This is the single parameter passed into each factory during creation.
    *
    * @return The backend configuration.
    */
  def configurationDescriptor: BackendConfigurationDescriptor

  /**
    * Returns the main engine for async execution.
    *
    * @return the main engine for async execution.
    */
  def asyncJobExecutionActorClass: Class[_ <: StandardAsyncExecutionActor]

  /**
    * Returns the key to use for storing and looking up the job id.
    *
    * @return the key to use for storing and looking up the job id.
    */
  def jobIdKey: String

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationDataOption: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef,
                                      backendSingletonActor: Option[ActorRef]): Props = {
    val params = DefaultStandardSyncExecutionActorParams(jobIdKey, serviceRegistryActor,
      jobDescriptor, configurationDescriptor, initializationDataOption, asyncJobExecutionActorClass)
    Props(new StandardSyncExecutionActor(params)).withDispatcher(Dispatcher.BackendDispatcher)
  }
}
