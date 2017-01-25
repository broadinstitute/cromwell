package cromwell.backend.standard

import java.nio.file.Path

import akka.actor.ActorRef
import cromwell.backend.callcaching.CacheHitDuplicating
import cromwell.backend.{BackendCacheHitCopyingActor, BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptor}

/**
  * Trait of parameters passed to a StandardCacheHitCopyingActor.
  */
trait StandardCacheHitCopyingActorParams {
  def jobDescriptor: BackendJobDescriptor

  def backendInitializationDataOption: Option[BackendInitializationData]

  def serviceRegistryActor: ActorRef

  def configurationDescriptor: BackendConfigurationDescriptor
}

/** A default implementation of the cache hit copying params. */
case class DefaultStandardCacheHitCopyingActorParams
(
  override val jobDescriptor: BackendJobDescriptor,
  override val backendInitializationDataOption: Option[BackendInitializationData],
  override val serviceRegistryActor: ActorRef,
  override val configurationDescriptor: BackendConfigurationDescriptor
) extends StandardCacheHitCopyingActorParams

/**
  * Standard implementation of a BackendCacheHitCopyingActor.
  */
trait StandardCacheHitCopyingActor extends BackendCacheHitCopyingActor with CacheHitDuplicating
  with StandardCachingActorHelper {

  def standardParams: StandardCacheHitCopyingActorParams

  override lazy val jobDescriptor: BackendJobDescriptor = standardParams.jobDescriptor
  override lazy val backendInitializationDataOption: Option[BackendInitializationData] =
    standardParams.backendInitializationDataOption
  override lazy val serviceRegistryActor: ActorRef = standardParams.serviceRegistryActor
  override lazy val configurationDescriptor: BackendConfigurationDescriptor = standardParams.configurationDescriptor
  override lazy val destinationCallRootPath: Path = jobPaths.callRoot
  override lazy val destinationJobDetritusPaths: Map[String, Path] = jobPaths.detritusPaths
}
