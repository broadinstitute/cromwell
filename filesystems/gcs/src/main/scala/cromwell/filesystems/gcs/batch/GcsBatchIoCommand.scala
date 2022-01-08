package cromwell.filesystems.gcs.batch

import cats.implicits.catsSyntaxValidatedId
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpHeaders
import com.google.api.services.storage.StorageRequest
import com.google.api.services.storage.model.{Objects, RewriteResponse, StorageObject}
import com.google.cloud.storage.BlobId
import common.util.StringUtil._
import common.validation.ErrorOr.ErrorOr
import cromwell.core.io._
import cromwell.filesystems.gcs._
import mouse.all._

import scala.collection.JavaConverters._
import scala.util.Try

/**
  * NOTE: the setUserProject commented out code disables uses of requester pays bucket
  * To re-enable, uncomment the code (and possibly adjust if necessary)
  *
  * Io commands with GCS paths and some logic enabling batching of request.
  * @tparam T Return type of the IoCommand
  * @tparam U Return type of the Google response
  */
sealed trait GcsBatchIoCommand[T, U] extends IoCommand[T] {
  /**
    * StorageRequest operation to be executed by this command
    */
  def operation: StorageRequest[U]

  /**
    * Maps the google response of type U to the Cromwell Io response of type T
    */
  protected def mapGoogleResponse(response: U): ErrorOr[T]

  /**
    * Method called in the success callback of a batched request to decide what to do next.
    * Returns an `ErrorOr[Either[T, GcsBatchIoCommand[T, U]]]`
    * Within the Either:
    *   Left(value) means the command is complete, and the result can be sent back to the sender.
    *   Right(newCommand) means the command is not complete and needs another request to be executed.
    * Most commands will reply with Left(value).
    */
  def onSuccess(response: U, httpHeaders: HttpHeaders): ErrorOr[Either[T, GcsBatchIoCommand[T, U]]] = {
    logIOMsgOverLimit(s"GcsBatchIoCommand.onSuccess '${response.toPrettyElidedString(limit = 1000)}'")
    mapGoogleResponse(response) map Left.apply
  }

  /**
    * Override to handle a failure differently and potentially return a successful response.
    */
  def onFailure(googleJsonError: GoogleJsonError, httpHeaders: HttpHeaders): Option[Either[T, GcsBatchIoCommand[T, U]]] = None

  /**
    * Use to signal that the request has failed because the user project was not set and that it can be retried with it.
    */
  def withUserProject: GcsBatchIoCommand[T, U]
}

sealed trait SingleFileGcsBatchIoCommand[T, U] extends GcsBatchIoCommand[T, U] with SingleFileIoCommand[T] {
  override def file: GcsPath
  //noinspection MutatorLikeMethodIsParameterless
  def setUserProject: Boolean
  def userProject: String = setUserProject.option(file.projectId).orNull
}

case class GcsBatchCopyCommand(
                                override val source: GcsPath,
                                sourceBlob: BlobId,
                                override val destination: GcsPath,
                                destinationBlob: BlobId,
                                rewriteToken: Option[String] = None,
                                setUserProject: Boolean = false
                              )
  extends IoCopyCommand(source, destination) with GcsBatchIoCommand[Unit, RewriteResponse] {
  override def commandDescription: String = s"GcsBatchCopyCommand source '$source' destination '$destination' " +
    s"setUserProject '$setUserProject' rewriteToken '$rewriteToken'"

  override def operation: StorageRequest[RewriteResponse] = {
    val rewriteOperation = source.apiStorage.objects()
      .rewrite(sourceBlob.getBucket, sourceBlob.getName, destinationBlob.getBucket, destinationBlob.getName, null)
      .setUserProject(setUserProject.option(source.projectId).orNull)

    // Set the rewrite token if present
    rewriteToken foreach rewriteOperation.setRewriteToken
    rewriteOperation
  }

  /**
    * Clone this command with the give rewrite token
    */
  def withRewriteToken(rewriteToken: String): GcsBatchCopyCommand = copy(rewriteToken = Option(rewriteToken))

  override def onSuccess(response: RewriteResponse, httpHeaders: HttpHeaders): ErrorOr[Either[Unit, GcsBatchCopyCommand]] = {
    if (response.getDone) {
      mapGoogleResponse(response) map Left.apply
    } else {
      Right(withRewriteToken(response.getRewriteToken)).validNel
    }
  }

  override def mapGoogleResponse(response: RewriteResponse): ErrorOr[Unit] = ().validNel

  override def withUserProject: GcsBatchCopyCommand = this.copy(setUserProject = true)
}

object GcsBatchCopyCommand {
  def forPaths(source: GcsPath, destination: GcsPath): Try[GcsBatchCopyCommand] = {
    for {
      sourceBlob <- source.objectBlobId
      destinationBlob <- destination.objectBlobId
    } yield GcsBatchCopyCommand(source, sourceBlob, destination, destinationBlob)
  }
}

case class GcsBatchDeleteCommand(
                                  override val file: GcsPath,
                                  blob: BlobId,
                                  override val swallowIOExceptions: Boolean,
                                  setUserProject: Boolean = false
                                ) extends IoDeleteCommand(file, swallowIOExceptions) with SingleFileGcsBatchIoCommand[Unit, Void] {
  override def operation: StorageRequest[Void] = {
    file.apiStorage.objects().delete(blob.getBucket, blob.getName).setUserProject(userProject)
  }

  override def mapGoogleResponse(response: Void): ErrorOr[Unit] = ().validNel

  override def onFailure(googleJsonError: GoogleJsonError,
                         httpHeaders: HttpHeaders,
                        ): Option[Either[Unit, GcsBatchDeleteCommand]] = {
    if (swallowIOExceptions) Option(Left(())) else None
  }
  override def withUserProject: GcsBatchDeleteCommand = this.copy(setUserProject = true)

  override def commandDescription: String = s"GcsBatchDeleteCommand file '$file' swallowIOExceptions " +
    s"'$swallowIOExceptions' setUserProject '$setUserProject'"
}

object GcsBatchDeleteCommand {
  def forPath(file: GcsPath, swallowIOExceptions: Boolean): Try[GcsBatchDeleteCommand] = {
    file.objectBlobId.map(GcsBatchDeleteCommand(file, _, swallowIOExceptions))
  }
}

/**
  * Base trait for commands that use the objects.get() operation. (e.g: size, crc32, ...)
  */
sealed trait GcsBatchGetCommand[T] extends SingleFileGcsBatchIoCommand[T, StorageObject] {
  def file: GcsPath
  def blob: BlobId
  override def operation: StorageRequest[StorageObject] = {
    file.apiStorage.objects().get(blob.getBucket, blob.getName).setUserProject(userProject)
  }
}

case class GcsBatchSizeCommand(override val file: GcsPath,
                               override val blob: BlobId,
                               setUserProject: Boolean = false,
                              ) extends IoSizeCommand(file) with GcsBatchGetCommand[Long] {
  override def mapGoogleResponse(response: StorageObject): ErrorOr[Long] = {
    Option(response.getSize) match {
      case None => s"'${file.pathAsString}' in project '${file.projectId}' returned null size".invalidNel
      case Some(size) => size.longValue().validNel
    }
  }

  override def withUserProject: GcsBatchSizeCommand = this.copy(setUserProject = true)

  override def commandDescription: String = s"GcsBatchSizeCommand file '$file' setUserProject '$setUserProject'"
}

object GcsBatchSizeCommand {
  def forPath(file: GcsPath): Try[GcsBatchSizeCommand] = {
    file.objectBlobId.map(GcsBatchSizeCommand(file, _))
  }
}

case class GcsBatchCrc32Command(override val file: GcsPath,
                                override val blob: BlobId,
                                setUserProject: Boolean = false,
                               ) extends IoHashCommand(file) with GcsBatchGetCommand[String] {
  override def mapGoogleResponse(response: StorageObject): ErrorOr[String] = {
    Option(response.getCrc32c) match {
      case None => s"'${file.pathAsString}' in project '${file.projectId}' returned null CRC32C checksum".invalidNel
      case Some(crc32c) => crc32c.validNel
    }
  }

  override def withUserProject: GcsBatchCrc32Command = this.copy(setUserProject = true)

  override def commandDescription: String = s"GcsBatchCrc32Command file '$file' setUserProject '$setUserProject'"
}

object GcsBatchCrc32Command {
  def forPath(file: GcsPath): Try[GcsBatchCrc32Command] = {
    file.objectBlobId.map(GcsBatchCrc32Command(file, _))
  }
}

case class GcsBatchTouchCommand(override val file: GcsPath,
                                override val blob: BlobId,
                                setUserProject: Boolean = false,
                               ) extends IoTouchCommand(file) with GcsBatchGetCommand[Unit] {
  override def mapGoogleResponse(response: StorageObject): ErrorOr[Unit] = ().validNel

  override def withUserProject: GcsBatchTouchCommand = this.copy(setUserProject = true)

  override def commandDescription: String = s"GcsBatchTouchCommand file '$file' setUserProject '$setUserProject'"
}

object GcsBatchTouchCommand {
  def forPath(file: GcsPath): Try[GcsBatchTouchCommand] = {
    file.objectBlobId.map(GcsBatchTouchCommand(file, _))
  }
}

/*
 * The only reliable way to know if a path represents a GCS "directory" is to list objects inside of it.
 * Specifically, list objects that have this path as a prefix. Since we don't really care about what's inside here,
 * set max results to 1 to avoid unnecessary payload.
 */
case class GcsBatchIsDirectoryCommand(override val file: GcsPath,
                                      blob: BlobId,
                                      setUserProject: Boolean = false,
                                     )
  extends IoIsDirectoryCommand(file) with SingleFileGcsBatchIoCommand[Boolean, Objects] {
  override def operation: StorageRequest[Objects] = {
    file.apiStorage.objects().list(blob.getBucket).setPrefix(blob.getName.ensureSlashed).setMaxResults(1L).setUserProject(userProject)
  }

  override def mapGoogleResponse(response: Objects): ErrorOr[Boolean] = {
    // Wrap in an Option because getItems can (always ?) return null if there are no objects
    Option(response.getItems).map(_.asScala).exists(_.nonEmpty).validNel
  }
  override def withUserProject: GcsBatchIsDirectoryCommand = this.copy(setUserProject = true)
  override def commandDescription: String = s"GcsBatchIsDirectoryCommand file '$file' setUserProject '$setUserProject'"
}

object GcsBatchIsDirectoryCommand {
  def forPath(file: GcsPath): Try[GcsBatchIsDirectoryCommand] = {
    file.bucketOrObjectBlobId.map(GcsBatchIsDirectoryCommand(file, _))
  }
}

case class GcsBatchExistsCommand(override val file: GcsPath,
                                 override val blob: BlobId,
                                 setUserProject: Boolean = false,
                                ) extends IoExistsCommand(file) with GcsBatchGetCommand[Boolean] {
  override def mapGoogleResponse(response: StorageObject): ErrorOr[Boolean] = true.validNel

  override def onFailure(googleJsonError: GoogleJsonError, httpHeaders: HttpHeaders): Option[Either[Boolean, GcsBatchIoCommand[Boolean, StorageObject]]] = {
    // If the object can't be found, don't fail the request but just return false as we were testing for existence
    if (googleJsonError.getCode == 404) Option(Left(false)) else None
  }
  override def withUserProject: GcsBatchExistsCommand = this.copy(setUserProject = true)
  override def commandDescription: String = s"GcsBatchExistsCommand file '$file' setUserProject '$setUserProject'"
}

object GcsBatchExistsCommand {
  def forPath(file: GcsPath): Try[GcsBatchExistsCommand] = {
    file.objectBlobId.map(GcsBatchExistsCommand(file, _))
  }
}

/** A GcsBatchIoCommand for use in tests. */
case object ExceptionSpewingGcsBatchIoCommand extends GcsBatchIoCommand[Unit, Void] {
  override lazy val name: String = "exception"

  override def operation: Nothing = throw new UnsupportedOperationException("operation is not supported")

  override def withUserProject: Nothing = throw new UnsupportedOperationException("withUserProject is not supported")

  override def mapGoogleResponse(response: Void): Nothing =
    sys.error("Ill behaved code that throws in mapGoogleResponse")

  override def onFailure(googleJsonError: GoogleJsonError, httpHeaders: HttpHeaders): Nothing =
    sys.error("Ill behaved code that throws in onFailure")

  override def commandDescription: String = s"ExceptionSpewingGcsBatchIoCommand (mock class for tests)"
}

/** A GcsBatchIoCommand for use in tests. */
case class ErrorReturningGcsBatchIoCommand() extends GcsBatchIoCommand[Unit, Unit] {
  override def operation: Nothing = throw new UnsupportedOperationException("operation is not supported")

  override def mapGoogleResponse(response: Unit): ErrorOr[Unit] =
    "Well behaved code that returns an error in mapGoogleResponse".invalidNel

  override def withUserProject: Nothing = throw new UnsupportedOperationException("withUserProject is not supported")

  /* Better to return a java.lang.Exception than try to shutdown the JVM with a java.lang.Error like this was before. */
  override def name: Nothing = throw new UnsupportedOperationException("a ErrorReturningGcsBatchIoCommand has no name")

  override def commandDescription: String = s"ErrorReturningGcsBatchIoCommand (mock class for tests)"
}
