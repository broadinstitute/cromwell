package cromwell.filesystems.drs

import java.io.ByteArrayInputStream
import java.nio.channels.ReadableByteChannel

import akka.actor.ActorSystem
import cats.data.Validated.{Invalid, Valid}
import cats.effect.IO
import cloud.nio.impl.drs.{DrsCloudNioFileSystemProvider, GcsFilePath, MarthaResponse}
import com.google.api.services.oauth2.Oauth2Scopes
import com.google.auth.oauth2.GoogleCredentials
import com.google.cloud.storage.Storage.BlobGetOption
import com.google.cloud.storage.{Blob, StorageException, StorageOptions}
import com.typesafe.config.Config
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.WorkflowOptions
import cromwell.core.path.{PathBuilder, PathBuilderFactory}
import org.apache.http.impl.client.HttpClientBuilder

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

  private lazy val httpClientBuilder = HttpClientBuilder.create()

  private val GcsScheme = "gs"


  private def gcsInputStream(gcsFile: GcsFilePath,
                             serviceAccountJson: String,
                             requesterPaysProjectIdOption: Option[String]): IO[ReadableByteChannel] = {
    val credentials = GoogleCredentials.fromStream(new ByteArrayInputStream(serviceAccountJson.getBytes()))
    val storage = StorageOptions.newBuilder().setCredentials(credentials).build().getService

    IO.delay {
      val blob = storage.get(gcsFile.bucket, gcsFile.file)
      blob.reader()
    } handleErrorWith { throwable =>
      (requesterPaysProjectIdOption, throwable) match {
        case (Some(requesterPaysProjectId), storageException: StorageException)
          if storageException.getMessage == "Bucket is requester pays bucket but no user project provided." =>
          IO.delay {
            val blob = storage.get(gcsFile.bucket, gcsFile.file, BlobGetOption.userProject(requesterPaysProjectId))
            blob.reader(Blob.BlobSourceOption.userProject(requesterPaysProjectId))
          }
        case _ => IO.raiseError(throwable)
      }
    }
  }


  private def inputReadChannel(url: String,
                               urlScheme: String,
                               serviceAccountJson: String,
                               requesterPaysProjectIdOption: Option[String]): IO[ReadableByteChannel] =  {
    urlScheme match {
      case GcsScheme =>
        val Array(bucket, fileToBeLocalized) = url.replace(s"$GcsScheme://", "").split("/", 2)
        gcsInputStream(GcsFilePath(bucket, fileToBeLocalized), serviceAccountJson, requesterPaysProjectIdOption)
      case otherScheme => IO.raiseError(new UnsupportedOperationException(s"DRS currently doesn't support reading files for $otherScheme."))
    }
  }


  private def drsReadInterpreter(options: WorkflowOptions, requesterPaysProjectIdOption: Option[String])
                                (marthaResponse: MarthaResponse): IO[ReadableByteChannel] = {
    val serviceAccountJsonIo = marthaResponse.googleServiceAccount match {
      case Some(googleSA) => IO.pure(googleSA.data.toString)
      case None => IO.fromEither(options.get(GoogleAuthMode.UserServiceAccountKey).toEither)
    }

    for {
      serviceAccountJson <- serviceAccountJsonIo
      //Currently, Martha only supports resolving DRS paths to GCS paths
      url <- IO.fromEither(DrsResolver.extractUrlForScheme(marthaResponse.dos.data_object.urls, GcsScheme))
      readableByteChannel <- inputReadChannel(url, GcsScheme, serviceAccountJson, requesterPaysProjectIdOption)
    } yield readableByteChannel
  }


  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[PathBuilder] = {
    val marthaScopes = List(
      // Profile and Email scopes are requirements for interacting with Martha v2
      Oauth2Scopes.USERINFO_EMAIL,
      Oauth2Scopes.USERINFO_PROFILE
    )
    val authCredentials = googleAuthMode.credentials(options.get(_).get, marthaScopes)

    // Unlike PAPI we're not going to fall back to a "default" project from the backend config.
    // ONLY use the project id from the User Service Account for requester pays
    val requesterPaysProjectIdOption = options.get("google_project").toOption

    Future(DrsPathBuilder(
      new DrsCloudNioFileSystemProvider(
        singletonConfig.config,
        authCredentials,
        httpClientBuilder,
        drsReadInterpreter(options, requesterPaysProjectIdOption)
      ),
      requesterPaysProjectIdOption,
    ))
  }
}



case class UrlNotFoundException(scheme: String) extends Exception(s"No $scheme url associated with given DRS path.")
