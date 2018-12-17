package cromwell.filesystems.drs

import akka.actor.ActorSystem
import cats.data.Validated.{Invalid, Valid}
import cloud.nio.impl.drs.DrsCloudNioFileSystemProvider
import com.google.api.services.oauth2.Oauth2Scopes
import com.typesafe.config.Config
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.WorkflowOptions
import cromwell.core.path.{PathBuilder, PathBuilderFactory}

import scala.concurrent.{ExecutionContext, Future}
import scala.collection.JavaConverters._


/**
  * Cromwell Wrapper around DrsFileSystems to load the configuration.
  * This class is used as the global configuration class in the drs filesystem
  */
class DrsFileSystemConfig(val config: Config)


class DrsPathBuilderFactory(globalConfig: Config, instanceConfig: Config, singletonConfig: DrsFileSystemConfig) extends PathBuilderFactory {

  private lazy val googleConfiguration: GoogleConfiguration = GoogleConfiguration(globalConfig)
  private lazy val scheme = instanceConfig.getString("auth")
  private lazy val googleAuthMode = googleConfiguration.auth(scheme) match {
    case Valid(auth) => auth
    case Invalid(error) => throw new RuntimeException(s"Error while instantiating DRS path builder factory. Errors: ${error.toString}")
  }


  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[PathBuilder] = {
    val marthaScopes = List(
      // Profile and Email scopes are requirements for interacting with Martha v2
      Oauth2Scopes.USERINFO_EMAIL,
      Oauth2Scopes.USERINFO_PROFILE
    ).asJavaCollection
    val authCredentials = googleAuthMode.credentials((key: String) => options.get(key).get, marthaScopes)

    Future.successful(DrsPathBuilder(new DrsCloudNioFileSystemProvider(singletonConfig.config, authCredentials)))
  }
}
