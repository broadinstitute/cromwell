package cromwell.filesystems.drs

import cats.effect.IO
import cloud.nio.impl.drs.DrsResolverResponseSupport.getGcsBucketAndName
import cloud.nio.impl.drs.{AccessUrl, DrsPathResolver, DrsResolverResponse, SADataObject}
import com.google.api.services.storage.StorageScopes
import com.google.auth.oauth2.OAuth2Credentials
import com.google.cloud.storage.Storage.BlobGetOption
import com.google.cloud.storage.{Blob, StorageException, StorageOptions}
import cromwell.cloudsupport.gcp.auth.{GoogleAuthMode, UserServiceAccountMode}
import cromwell.core.WorkflowOptions

import java.nio.channels.ReadableByteChannel

trait DrsReader {
  def read(): IO[ReadableByteChannel]
}

object DrsReader {
  def reader(googleAuthMode: Option[GoogleAuthMode],
             options: WorkflowOptions,
             requesterPaysProjectIdOption: Option[String],
             drsPathResolver: DrsPathResolver,
             drsResolverResponse: DrsResolverResponse): IO[DrsReader] = {
    (drsResolverResponse.accessUrl, drsResolverResponse.gsUri, googleAuthMode) match {
      case (Some(accessUrl), _, _) =>
        IO.pure(AccessUrlReader(drsPathResolver, accessUrl))
      case (_, Some(gcsPath), Some(authMode)) =>
        IO.pure(GcsReader(
          authMode,
          options,
          requesterPaysProjectIdOption,
          gcsPath,
          drsResolverResponse.googleServiceAccount,
        ))
      case (_, Some(_), _) =>
        IO.raiseError(new RuntimeException("GCS URI found in the DRS Resolver response, but no Google auth found!"))
      case _ =>
        IO.raiseError(new RuntimeException(DrsPathResolver.ExtractUriErrorMsg))
    }
  }

  def readInterpreter(googleAuthMode: Option[GoogleAuthMode],
                      options: WorkflowOptions,
                      requesterPaysProjectIdOption: Option[String])
                     (drsPathResolver: DrsPathResolver,
                      drsResolverResponse: DrsResolverResponse): IO[ReadableByteChannel] = {
    for {
      reader <- reader(googleAuthMode, options, requesterPaysProjectIdOption, drsPathResolver, drsResolverResponse)
      channel <- reader.read()
    } yield channel
  }
}

case class AccessUrlReader(drsPathResolver: DrsPathResolver, accessUrl: AccessUrl) extends DrsReader {
  override def read(): IO[ReadableByteChannel] = {
    drsPathResolver.openChannel(accessUrl)
  }
}

case class GcsReader(googleAuthMode: GoogleAuthMode,
                     options: WorkflowOptions,
                     requesterPaysProjectIdOption: Option[String],
                     gsUri: String,
                     googleServiceAccount: Option[SADataObject],
                    ) extends DrsReader {
  override def read(): IO[ReadableByteChannel] = {
    val readScopes = List(StorageScopes.DEVSTORAGE_READ_ONLY)
    val credentialsIo = googleServiceAccount match {
      case Some(googleSA) =>
        IO(
          UserServiceAccountMode("drs_resolver_service_account").credentials(
            Map(GoogleAuthMode.UserServiceAccountKey -> googleSA.data.noSpaces),
            readScopes,
          )
        )
      case None =>
        IO(googleAuthMode.credentials(options.get(_).get, readScopes))
    }

    for {
      credentials <- credentialsIo
      readableByteChannel <- gcsInputStream(gsUri, credentials, requesterPaysProjectIdOption)
    } yield readableByteChannel
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
              if storageException.getMessage.contains("requester pays bucket but no user project") =>
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
}
