package cromwell.backend.standard

import akka.actor.{Actor, ActorRef}
import cromwell.backend._
import cromwell.backend.callcaching.JobCachingActorHelper
import cromwell.backend.io.{JobPaths, WorkflowPaths}
import cromwell.backend.validation.{RuntimeAttributesValidation, ValidatedRuntimeAttributes}
import cromwell.core.logging.JobLogging
import cromwell.core.path.Path
import wdl4s.TaskCall

import scala.util.Try

/**
  * Extends the JobCachingActorHelper with standard implementations.
  *
  * Like the JobCachingActorHelper, this trait should be extended by a backend, and then __that__ trait should be mixed
  * into a async job execution actor, and cache hit copying actor.
  */
trait StandardCachingActorHelper extends JobCachingActorHelper {
  this: Actor with JobLogging =>

  def backendInitializationDataOption: Option[BackendInitializationData]

  /** Typed backend initialization. */
  def backendInitializationDataAs[A <: BackendInitializationData]: A =
    BackendInitializationData.as[A](backendInitializationDataOption)

  /**
    * Returns the service registry actor. Both the StandardAsyncExecutorActor and StandardCacheHitCopyingActor traits
    * implement this method.
    *
    * @return Paths to the job.
    */
  def serviceRegistryActor: ActorRef

  // So... JobPaths doesn't extend WorkflowPaths, but does contain a self-type
  lazy val workflowPaths: WorkflowPaths = jobPaths.asInstanceOf[WorkflowPaths]

  def getPath(str: String): Try[Path] = workflowPaths.getPath(str)

  /**
    * The workflow descriptor for this job. NOTE: This may be different than the workflow descriptor created in the
    * workflow initialization data. For example, sub workflows use a different workflow descriptor.
    */
  lazy val workflowDescriptor: BackendWorkflowDescriptor = jobDescriptor.workflowDescriptor

  lazy val call: TaskCall = jobDescriptor.key.call

  lazy val standardInitializationData: StandardInitializationData = BackendInitializationData.
    as[StandardInitializationData](backendInitializationDataOption)

  lazy val validatedRuntimeAttributes: ValidatedRuntimeAttributes = {
    val builder = standardInitializationData.runtimeAttributesBuilder
    builder.build(jobDescriptor.runtimeAttributes, jobLogger)
  }

  /**
    * Returns the paths to the job.
    *
    * @return Paths to the job.
    */
  lazy val jobPaths: JobPaths = standardInitializationData.workflowPaths.toJobPaths(jobDescriptor)

  /**
    * Returns the metadata key values to store before executing a job.
    *
    * @return the metadata key values to store before executing a job.
    */
  def startMetadataKeyValues: Map[String, Any] = {
    val runtimeAttributesMetadata = RuntimeAttributesValidation.toMetadataStrings(validatedRuntimeAttributes) map {
      case (key, value) => (s"${CallMetadataKeys.RuntimeAttributes}:$key", value)
    }

    val fileMetadata = jobPaths.metadataPaths

    runtimeAttributesMetadata ++ fileMetadata ++ nonStandardMetadata
  }

  /**
    * Returns any custom medatata for the backend.
    *
    * @return any custom metadata for the backend.
    */
  protected def nonStandardMetadata: Map[String, Any] = Map.empty
}
