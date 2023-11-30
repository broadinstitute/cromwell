package cromwell.filesystems.drs

import akka.actor.ActorSystem
import cats.data.Validated.{Invalid, Valid}
import cloud.nio.impl.drs.{AzureDrsCredentials, DrsCloudNioFileSystemProvider, GoogleOauthDrsCredentials}
import com.google.api.services.oauth2.Oauth2Scopes
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

class DrsPathBuilderFactory(globalConfig: Config, instanceConfig: Config, singletonConfig: DrsFileSystemConfig)
    extends PathBuilderFactory {

  private lazy val googleConfiguration: GoogleConfiguration = GoogleConfiguration(globalConfig)
  private lazy val scheme = instanceConfig.getString("auth")

  // For Azure support - this should be the UAMI client id
  private val dataAccessIdentityKey = "data_access_identity"

  override def withOptions(
    options: WorkflowOptions
  )(implicit as: ActorSystem, ec: ExecutionContext): Future[PathBuilder] =
    Future {
      val drsResolverScopes = List(
        // Profile and Email scopes are requirements for interacting with DRS Resolvers
        Oauth2Scopes.USERINFO_EMAIL,
        Oauth2Scopes.USERINFO_PROFILE
      )

      val (googleAuthMode, drsCredentials) = scheme match {
        case "azure" => (None, AzureDrsCredentials(options.get(dataAccessIdentityKey).toOption))
        case googleAuthScheme =>
          googleConfiguration.auth(googleAuthScheme) match {
            case Valid(auth) =>
              (
                Option(auth),
                GoogleOauthDrsCredentials(auth.credentials(options.get(_).get, drsResolverScopes),
                                          singletonConfig.config
                )
              )
            case Invalid(error) =>
              throw new RuntimeException(
                s"Error while instantiating DRS path builder factory. Errors: ${error.toString}"
              )
          }
      }

      // Unlike PAPI we're not going to fall back to a "default" project from the backend config.
      // ONLY use the project id from the User Service Account for requester pays
      val requesterPaysProjectIdOption = options.get("google_project").toOption

      /*
      `override_preresolve_for_test` is a workflow option to override the default `resolver.preresolve` specified in the
      global config. This is only used for testing purposes.
       */
      val preResolve: Boolean =
        options
          .getBoolean("override_preresolve_for_test")
          .toOption
          .getOrElse(
            singletonConfig.config
              .getBoolean("resolver.preresolve")
          )

      DrsPathBuilder(
        new DrsCloudNioFileSystemProvider(
          singletonConfig.config,
          drsCredentials,
          DrsReader.readInterpreter(googleAuthMode, options, requesterPaysProjectIdOption)
        ),
        requesterPaysProjectIdOption,
        preResolve
      )
    }
}

case class UrlNotFoundException(scheme: String) extends Exception(s"No $scheme url associated with given DRS path.")

case class DrsResolverResponseMissingKeyException(missingKey: String)
    extends Exception(s"The response from the DRS Resolver doesn't contain the key '$missingKey'.")
