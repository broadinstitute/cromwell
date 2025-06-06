package cromwell.backend.standard

import akka.actor.ActorRef
import cromwell.backend.{
  BackendConfigurationDescriptor,
  BackendInitializationData,
  BackendJobDescriptor,
  MinimumRuntimeSettings
}

/**
  * Base trait for params passed to both the sync and async backend actors.
  */
trait StandardJobExecutionActorParams {

  /** The service registry actor for key/value and metadata. */
  def serviceRegistryActor: ActorRef

  /** Actor able to handle IO requests asynchronously */
  def ioActor: ActorRef

  /** The descriptor of this job. */
  def jobDescriptor: BackendJobDescriptor

  /** The global and backend configuration. */
  def configurationDescriptor: BackendConfigurationDescriptor

  /** Any backend initialization data. */
  def backendInitializationDataOption: Option[BackendInitializationData]

  /** The key for this job. */
  def jobIdKey: String

  /** The singleton actor. */
  def backendSingletonActorOption: Option[ActorRef]

  /** Singleton actor for recording when hog group runs into quota exhaustion */
  def groupMetricsActor: ActorRef

  /** The default settings for runtime Environment passed to CWL expressions when not specified in the Resource Requirements */
  val minimumRuntimeSettings: MinimumRuntimeSettings
}
