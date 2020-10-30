package cromwell.filesystems.gcs.batch

import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpHeaders
import com.google.api.services.storage.StorageRequest
import com.google.api.services.storage.model.{Objects, RewriteResponse, StorageObject}
import common.util.StringUtil._
import cromwell.core.io._
import cromwell.filesystems.gcs._
import mouse.all._

import scala.collection.JavaConverters._

/**
  * NOTE: the setUserProject commented out code disables uses of requester pays bucket
  * To re-enable, uncomment the code (and possibly adjust if necessary)
  *
  * Io commands with GCS paths and some logic enabling batching of request.
  * @tparam T Return type of the IoCommand
  * @tparam U Return type of the Google response
  */
trait GcsBatchIoCommand[T, U] extends IoCommand[T] {
  /**
    * StorageRequest operation to be executed by this command
    */
  def operation: StorageRequest[U]

  /**
    * Maps the google response of type U to the Cromwell Io response of type T
    */
  protected def mapGoogleResponse(response: U): T

  /**
    * Method called in the success callback of a batched request to decide what to do next.
    * Returns an Either[T, GcsBatchIoCommand[T, U]]
    *   Left(value) means the command is complete, and the result can be sent back to the sender.
    *   Right(newCommand) means the command is not complete and needs another request to be executed.
    * Most commands will reply with Left(value).
    */
  def onSuccess(response: U, httpHeaders: HttpHeaders): Either[T, GcsBatchIoCommand[T, U]] = {
    logIOMsgOverLimit(s"GcsBatchIoCommand.onSuccess '${response.toPrettyElidedString(limit = 1000)}'")
    Left(mapGoogleResponse(response))
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
  def setUserProject: Boolean
  def userProject = setUserProject.option(file.projectId).orNull
}

case class GcsBatchCopyCommand(
                                override val source: GcsPath,
                                override val destination: GcsPath,
                                override val overwrite: Boolean,
                                rewriteToken: Option[String] = None,
                                setUserProject: Boolean = false
                              ) extends IoCopyCommand(source, destination, overwrite) with GcsBatchIoCommand[Unit, RewriteResponse] {
  val sourceBlob = source.blob
  val destinationBlob = destination.blob

  override def commandDescription: String = s"GcsBatchCopyCommand source '$source' destination '$destination' " +
    s"overwrite '$overwrite' setUserProject '$setUserProject'"

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
  def withRewriteToken(rewriteToken: String) = copy(rewriteToken = Option(rewriteToken))

  override def onSuccess(response: RewriteResponse, httpHeaders: HttpHeaders) = {
    if (response.getDone) super.onSuccess(response, httpHeaders)
    else {
      Right(withRewriteToken(response.getRewriteToken))
    }
  }

  override def mapGoogleResponse(response: RewriteResponse): Unit = ()
  override def withUserProject = this.copy(setUserProject = true)
}

case class GcsBatchDeleteCommand(
                                  override val file: GcsPath,
                                  override val swallowIOExceptions: Boolean,
                                  setUserProject: Boolean = false
                                ) extends IoDeleteCommand(file, swallowIOExceptions) with SingleFileGcsBatchIoCommand[Unit, Void] {
  private val blob = file.blob
  def operation = {
    file.apiStorage.objects().delete(blob.getBucket, blob.getName).setUserProject(userProject)
  }
  override protected def mapGoogleResponse(response: Void): Unit = ()
  override def onFailure(googleJsonError: GoogleJsonError, httpHeaders: HttpHeaders) = {
    if (swallowIOExceptions) Option(Left(())) else None
  }
  override def withUserProject = this.copy(setUserProject = true)

  override def commandDescription: String = s"GcsBatchDeleteCommand file '$file' swallowIOExceptions " +
    s"'$swallowIOExceptions' setUserProject '$setUserProject'"
}

/**
  * Base trait for commands that use the objects.get() operation. (e.g: size, crc32, ...)
  */
sealed trait GcsBatchGetCommand[T] extends SingleFileGcsBatchIoCommand[T, StorageObject] {
  def file: GcsPath
  private val blob = file.blob
  override def operation: StorageRequest[StorageObject] = {
    file.apiStorage.objects().get(blob.getBucket, blob.getName).setUserProject(userProject)
  }
}

case class GcsBatchSizeCommand(override val file: GcsPath, setUserProject: Boolean = false) extends IoSizeCommand(file) with GcsBatchGetCommand[Long] {
  override def mapGoogleResponse(response: StorageObject): Long = response.getSize.longValue()
  override def withUserProject = this.copy(setUserProject = true)
  override def commandDescription: String = s"GcsBatchSizeCommand file '$file' setUserProject '$setUserProject'"
}

case class GcsBatchCrc32Command(override val file: GcsPath, setUserProject: Boolean = false) extends IoHashCommand(file) with GcsBatchGetCommand[String] {
  override def mapGoogleResponse(response: StorageObject): String = response.getCrc32c
  override def withUserProject = this.copy(setUserProject = true)
  override def commandDescription: String = s"GcsBatchCrc32Command file '$file' setUserProject '$setUserProject'"
}

case class GcsBatchTouchCommand(override val file: GcsPath, setUserProject: Boolean = false) extends IoTouchCommand(file) with GcsBatchGetCommand[Unit] {
  override def mapGoogleResponse(response: StorageObject): Unit = ()
  override def withUserProject = this.copy(setUserProject = true)
  override def commandDescription: String = s"GcsBatchTouchCommand file '$file' setUserProject '$setUserProject'"
}

/*
 * The only reliable way to know if a path represents a GCS "directory" is to list objects inside of it.
 * Specifically, list objects that have this path as a prefix. Since we don't really care about what's inside here,
 * set max results to 1 to avoid unnecessary payload.
 */
case class GcsBatchIsDirectoryCommand(override val file: GcsPath, setUserProject: Boolean = false) extends IoIsDirectoryCommand(file) with SingleFileGcsBatchIoCommand[Boolean, Objects] {
  private val blob = file.blob
  override def operation: StorageRequest[Objects] = {
    file.apiStorage.objects().list(blob.getBucket).setPrefix(blob.getName.ensureSlashed).setMaxResults(1L).setUserProject(userProject)
  }
  override def mapGoogleResponse(response: Objects): Boolean = {
    // Wrap in an Option because getItems can (always ?) return null if there are no objects
    Option(response.getItems).map(_.asScala).exists(_.nonEmpty)
  }
  override def withUserProject = this.copy(setUserProject = true)
  override def commandDescription: String = s"GcsBatchIsDirectoryCommand file '$file' setUserProject '$setUserProject'"
}

case class GcsBatchExistsCommand(override val file: GcsPath, setUserProject: Boolean = false) extends IoExistsCommand(file) with GcsBatchGetCommand[Boolean] {
  override def mapGoogleResponse(response: StorageObject): Boolean = true
  override def onFailure(googleJsonError: GoogleJsonError, httpHeaders: HttpHeaders) = {
    // If the object can't be found, don't fail the request but just return false as we were testing for existence
    if (googleJsonError.getCode == 404) Option(Left(false)) else None
  }
  override def withUserProject = this.copy(setUserProject = true)
  override def commandDescription: String = s"GcsBatchExistsCommand file '$file' setUserProject '$setUserProject'"
}
