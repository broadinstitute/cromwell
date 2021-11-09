package cromwell.filesystems.drs

import akka.actor.ActorSystem
import cats.data.Validated.{Invalid, Valid}
import cloud.nio.impl.drs.DrsCloudNioFileSystemProvider
import com.google.api.services.oauth2.Oauth2Scopes
import com.google.cloud.NoCredentials
import com.typesafe.config.Config
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.WorkflowOptions
import cromwell.core.path.{PathBuilder, PathBuilderFactory}

import scala.concurrent.{ExecutionContext, Future}

/**
  * Cromwell Wrapper around DrsFileSystems to load the configuration.
  * This class is used as the global configuration class in the drs filesystem
  */
class DrsFileSystemConfig(val config: Config)


class DrsPathBuilderFactory(globalConfig: Config, instanceConfig: Config, singletonConfig: DrsFileSystemConfig) extends PathBuilderFactory {

  private lazy val googleConfiguration: GoogleConfiguration = GoogleConfiguration(globalConfig)
  private lazy val scheme = instanceConfig.getString("auth")

  // TODO
  private lazy val googleAuthMode = googleConfiguration.auth(scheme) match {
    case Valid(auth) => Some(auth)
    case Invalid(_) => None//throw new RuntimeException(s"Error while instantiating DRS path builder factory. Errors: ${error.toString}")
  }

  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[PathBuilder] = {
    Future {
      val marthaScopes = List(
        // Profile and Email scopes are requirements for interacting with Martha
        Oauth2Scopes.USERINFO_EMAIL,
        Oauth2Scopes.USERINFO_PROFILE
      )
      val authCredentials =
        googleAuthMode
          .map(a => a.credentials(options.get(_).get, marthaScopes))
          .getOrElse(NoCredentials.getInstance())

      // Unlike PAPI we're not going to fall back to a "default" project from the backend config.
      // ONLY use the project id from the User Service Account for requester pays
      val requesterPaysProjectIdOption = options.get("google_project").toOption

      /*
      `override_preresolve_for_test` is a workflow option to override the default `martha.preresolve` specified in the
      global config. This is only used for testing purposes.
       */
      val preResolve: Boolean =
        options
          .getBoolean("override_preresolve_for_test")
          .toOption
          .getOrElse(
            singletonConfig
              .config
              .getBoolean("martha.preresolve")
          )

      DrsPathBuilder(
        new DrsCloudNioFileSystemProvider(
          singletonConfig.config,
          authCredentials,
          DrsReader.readInterpreter(googleAuthMode, options, requesterPaysProjectIdOption),
        ),
        requesterPaysProjectIdOption,
        preResolve,
      )
    }
  }
}

case class UrlNotFoundException(scheme: String) extends Exception(s"No $scheme url associated with given DRS path.")

case class MarthaResponseMissingKeyException(missingKey: String) extends Exception(s"The response from Martha doesn't contain the key '$missingKey'.")
