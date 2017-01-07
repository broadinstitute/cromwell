package cromwell.backend.standard

import akka.actor.ActorRef
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptor}

import scala.concurrent.Promise
import scala.language.existentials

/**
  * Base trait for params passed to both the sync and async backend actors.
  */
trait StandardJobExecutionActorParams {
  /** The service registry actor for key/value and metadata. */
  def serviceRegistryActor: ActorRef

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
}

/**
  * Extended trait for params passed to synchronous backend actors.
  */
trait StandardSyncExecutionActorParams extends StandardJobExecutionActorParams {
  /**
    * The class for creating an async backend.
    *
    * @see [[StandardSyncExecutionActor]]
    */
  def asyncJobExecutionActorClass: Class[_ <: StandardAsyncExecutionActor]
}

/** A default implementation of the sync params. */
case class DefaultStandardSyncExecutionActorParams
(
  override val jobIdKey: String,
  override val serviceRegistryActor: ActorRef,
  override val jobDescriptor: BackendJobDescriptor,
  override val configurationDescriptor: BackendConfigurationDescriptor,
  override val backendInitializationDataOption: Option[BackendInitializationData],
  override val backendSingletonActorOption: Option[ActorRef],
  override val asyncJobExecutionActorClass: Class[_ <: StandardAsyncExecutionActor]
) extends StandardSyncExecutionActorParams

/**
  * Extended trait for params passed to asynchronous backend actors.
  */
trait StandardAsyncExecutionActorParams extends StandardJobExecutionActorParams {
  /** The promise that will be completed when the async run is complete. */
  def completionPromise: Promise[BackendJobExecutionResponse]
}

/** A default implementation of the async params. */
case class DefaultStandardAsyncExecutionActorParams
(
  override val jobIdKey: String,
  override val serviceRegistryActor: ActorRef,
  override val jobDescriptor: BackendJobDescriptor,
  override val configurationDescriptor: BackendConfigurationDescriptor,
  override val backendInitializationDataOption: Option[BackendInitializationData],
  override val backendSingletonActorOption: Option[ActorRef],
  override val completionPromise: Promise[BackendJobExecutionResponse]
) extends StandardAsyncExecutionActorParams
