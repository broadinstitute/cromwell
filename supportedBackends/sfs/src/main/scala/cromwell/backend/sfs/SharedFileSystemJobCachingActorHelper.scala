package cromwell.backend.sfs

import akka.actor.{Actor, ActorRef}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.BackendInitializationData
import cromwell.backend.callcaching.JobCachingActorHelper
import cromwell.backend.io.{JobPathsWithDocker, WorkflowPathsBackendInitializationData}
import cromwell.backend.validation.{RuntimeAttributesValidation, ValidatedRuntimeAttributes}
import cromwell.core.logging.JobLogging
import cromwell.core.path.PathBuilder
import net.ceedubs.ficus.Ficus._

trait SharedFileSystemJobCachingActorHelper extends JobCachingActorHelper {
  this: Actor with JobLogging =>

  def backendInitializationDataOption: Option[BackendInitializationData]

  def serviceRegistryActor: ActorRef

  lazy val jobPaths =
    new JobPathsWithDocker(jobDescriptor.key, jobDescriptor.workflowDescriptor, configurationDescriptor.backendConfig)

  lazy val initializationData: SharedFileSystemBackendInitializationData = BackendInitializationData.
    as[SharedFileSystemBackendInitializationData](backendInitializationDataOption)

  lazy val validatedRuntimeAttributes: ValidatedRuntimeAttributes = {
    val builder = initializationData.runtimeAttributesBuilder
    builder.build(jobDescriptor.runtimeAttributes, jobLogger)
  }

  def startMetadataKeyValues: Map[String, Any] = {
    val runtimeAttributesMetadata = RuntimeAttributesValidation.extract(validatedRuntimeAttributes) map {
      case (key, value) => (s"runtimeAttributes:$key", value)
    }
    val fileMetadata = jobPaths.metadataPaths
    val otherMetadata = Map("cache:allowResultReuse" -> true)
    runtimeAttributesMetadata ++ fileMetadata ++ otherMetadata
  }

  lazy val sharedFileSystem = new SharedFileSystem {
    override val pathBuilders: List[PathBuilder] = {
      WorkflowPathsBackendInitializationData.pathBuilders(backendInitializationDataOption)
    }
    override lazy val sharedFileSystemConfig: Config = {
      configurationDescriptor.backendConfig.as[Option[Config]]("filesystems.local").getOrElse(ConfigFactory.empty())
    }
  }
}
