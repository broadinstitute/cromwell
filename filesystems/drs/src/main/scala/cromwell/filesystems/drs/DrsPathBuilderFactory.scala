package cromwell.filesystems.drs

import java.nio.channels.ReadableByteChannel

import akka.actor.ActorSystem
import cats.data.Validated.{Invalid, Valid}
import cats.effect.IO
import cloud.nio.impl.drs.MarthaResponseSupport._
import cloud.nio.impl.drs.{DrsCloudNioFileSystemProvider, SADataObject}
import com.google.api.services.oauth2.Oauth2Scopes
import com.google.api.services.storage.StorageScopes
import com.google.auth.oauth2.OAuth2Credentials
import com.google.cloud.storage.Storage.BlobGetOption
import com.google.cloud.storage.{Blob, StorageException, StorageOptions}
import com.typesafe.config.Config
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.{GoogleAuthMode, UserServiceAccountMode}
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
  private lazy val googleAuthMode = googleConfiguration.auth(scheme) match {
    case Valid(auth) => auth
    case Invalid(error) => throw new RuntimeException(s"Error while instantiating DRS path builder factory. Errors: ${error.toString}")
  }

  private def gcsInputStream(gcsFile: String,
                             credentials: OAuth2Credentials,
                             requesterPaysProjectIdOption: Option[String],
                            ): IO[ReadableByteChannel] = {
    for {
      storage <- IO(StorageOptions.newBuilder().setCredentials(credentials).build().getService)
      gcsBucketAndName <- IO(getGcsBucketAndName(gcsFile))
      (bucketName, objectName) = gcsBucketAndName
      readChannel <- IO(storage.get(bucketName, objectName).reader()) handleErrorWith {
        throwable =>
          (requesterPaysProjectIdOption, throwable) match {
            case (Some(requesterPaysProjectId), storageException: StorageException)
              if storageException.getMessage == "Bucket is requester pays bucket but no user project provided." =>
              IO(
                storage
                  .get(bucketName, objectName, BlobGetOption.userProject(requesterPaysProjectId))
                  .reader(Blob.BlobSourceOption.userProject(requesterPaysProjectId))
              )
            case _ => IO.raiseError(throwable)
          }
      }
    } yield readChannel
  }


  private def drsReadInterpreter(options: WorkflowOptions, requesterPaysProjectIdOption: Option[String])
                                (gsUri: Option[String], googleServiceAccount: Option[SADataObject])
  : IO[ReadableByteChannel] = {
    val readScopes = List(StorageScopes.DEVSTORAGE_READ_ONLY)
    val credentialsIo = googleServiceAccount match {
      case Some(googleSA) =>
        IO(
          UserServiceAccountMode("martha_service_account").credentials(
            Map(GoogleAuthMode.UserServiceAccountKey -> googleSA.data.noSpaces),
            readScopes,
          )
        )
      case None =>
        IO(googleAuthMode.credentials(options.get(_).get, readScopes))
    }

    for {
      credentials <- credentialsIo
      //Currently, Martha only supports resolving DRS paths to GCS paths
      gsUri <- IO.fromEither(gsUri.toRight(MarthaResponseMissingKeyException("gsUri")))
      readableByteChannel <- gcsInputStream(gsUri, credentials, requesterPaysProjectIdOption)
    } yield readableByteChannel
  }


  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[PathBuilder] = {
    val marthaScopes = List(
      // Profile and Email scopes are requirements for interacting with Martha
      Oauth2Scopes.USERINFO_EMAIL,
      Oauth2Scopes.USERINFO_PROFILE
    )
    val authCredentials = googleAuthMode.credentials(options.get(_).get, marthaScopes)

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

    Future(DrsPathBuilder(
      new DrsCloudNioFileSystemProvider(
        singletonConfig.config,
        authCredentials,
        drsReadInterpreter(options, requesterPaysProjectIdOption),
      ),
      requesterPaysProjectIdOption,
      preResolve,
    ))
  }
}



case class UrlNotFoundException(scheme: String) extends Exception(s"No $scheme url associated with given DRS path.")

case class MarthaResponseMissingKeyException(missingKey: String) extends Exception(s"The response from Martha doesn't contain the key '$missingKey'.")
