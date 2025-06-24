package cromwell.backend.impl.sfs.config

import com.typesafe.config.Config
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.backend.impl.sfs.config.ConfigConstants._
import cromwell.backend.sfs._
import cromwell.backend.standard.callcaching.StandardFileHashingActor
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.{BackendDockerConfiguration, DockerCredentials}
import net.ceedubs.ficus.Ficus._
import org.slf4j.{Logger, LoggerFactory}

import scala.util.{Success, Try}

/**
  * Builds a backend by reading the job control from the config.
  *
  * @param configurationDescriptor The config information.
  */
class ConfigBackendLifecycleActorFactory(val name: String, val configurationDescriptor: BackendConfigurationDescriptor)
    extends SharedFileSystemBackendLifecycleActorFactory {

  lazy val logger: Logger = LoggerFactory.getLogger(getClass)
  lazy val hashingStrategy: ConfigHashingStrategy =
    configurationDescriptor.backendConfig.as[Option[Config]](
      "filesystems.local.caching"
    ) map ConfigHashingStrategy.apply getOrElse ConfigHashingStrategy.defaultStrategy

  override lazy val initializationActorClass: Class[ConfigInitializationActor] = classOf[ConfigInitializationActor]

  override lazy val asyncExecutionActorClass: Class[_ <: ConfigAsyncJobExecutionActor] = {
    val runInBackground =
      configurationDescriptor.backendConfig.as[Option[Boolean]](RunInBackgroundConfig).getOrElse(false)
    if (runInBackground)
      classOf[BackgroundConfigAsyncJobExecutionActor]
    else
      classOf[DispatchedConfigAsyncJobExecutionActor]
  }

  override lazy val fileHashingActorClassOption: Option[Class[_ <: StandardFileHashingActor]] = Option(
    classOf[ConfigBackendFileHashingActor]
  )

  override def dockerHashCredentials(workflowDescriptor: BackendWorkflowDescriptor,
                                     initializationData: Option[BackendInitializationData]
  ): List[Any] =
    /*
    Heavily adapted from:
    https://github.com/broadinstitute/cromwell/blob/78/supportedBackends/google/pipelines/common/src/main/scala/cromwell/backend/google/pipelines/common/PipelinesApiBackendLifecycleActorFactory.scala#L71-L85

    Could also be moved into a "standard" location.
     */
    Try(BackendInitializationData.as[ConfigInitializationData](initializationData)) match {
      case Success(configInitializationData) =>
        val tokenFromWorkflowOptions =
          workflowDescriptor.workflowOptions
            .get(GoogleAuthMode.DockerCredentialsTokenKey)
            .toOption
        val effectiveToken =
          tokenFromWorkflowOptions
            .orElse(
              BackendDockerConfiguration
                .build(configurationDescriptor.backendConfig)
                .dockerCredentials
                .map(_.token)
            )

        val dockerCredentials: Option[DockerCredentials] = effectiveToken map { token =>
          // These credentials are being returned for hashing and all that matters in this context is the token
          // so just `None` the auth and key.
          val baseDockerCredentials = new DockerCredentials(token = token, authName = None, keyName = None)
          baseDockerCredentials
        }
        val googleCredentials = configInitializationData.googleRegistryCredentialsOption
        List(dockerCredentials, googleCredentials).flatten
      case _ => List.empty[Any]
    }
}
