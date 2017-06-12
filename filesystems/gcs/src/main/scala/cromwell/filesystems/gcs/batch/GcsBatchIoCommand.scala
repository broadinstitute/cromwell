package cromwell.filesystems.gcs.batch

import com.google.api.client.http.HttpHeaders
import com.google.api.services.storage.StorageRequest
import com.google.api.services.storage.model.{RewriteResponse, StorageObject}
import cromwell.core.io._
import cromwell.filesystems.gcs._

/**
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
  protected def mapGoogleResponse(response: U): T

  /**
    * Method called in the success callback of a batched request to decide what to do next.
    * Returns an Either[T, GcsBatchIoCommand[T, U]]
    *   Left(value) means the command is complete, and the result can be sent back to the sender.
    *   Right(newCommand) means the command is not complete and needs another request to be executed.
    * Most commands will reply with Left(value).
    */
  def onSuccess(response: U, httpHeaders: HttpHeaders): Either[T, GcsBatchIoCommand[T, U]] = {
    Left(mapGoogleResponse(response))
  }
}

case class GcsBatchCopyCommand(
                           override val source: GcsPath,
                           override val destination: GcsPath,
                           override val overwrite: Boolean,
                           rewriteToken: Option[String] = None
                         ) extends IoCopyCommand(source, destination, overwrite) with GcsBatchIoCommand[Unit, RewriteResponse] {
  val sourceBlob = source.blob
  val destinationBlob = destination.blob
  
  override def operation: StorageRequest[RewriteResponse] = {
    val rewriteOperation = source.apiStorage.objects().rewrite(sourceBlob.getBucket, sourceBlob.getName, destinationBlob.getBucket, destinationBlob.getName, null)
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
}

case class GcsBatchDeleteCommand(
                                  override val file: GcsPath,
                                  override val swallowIOExceptions: Boolean
                                ) extends IoDeleteCommand(file, swallowIOExceptions) with GcsBatchIoCommand[Unit, Void] {
  private val blob = file.blob
  def operation = file.apiStorage.objects().delete(blob.getBucket, blob.getName)
  override protected def mapGoogleResponse(response: Void): Unit = ()
}

/**
  * Base trait for commands that use the objects.get() operation. (e.g: size, crc32, ...)
  */
sealed trait GcsBatchGetCommand[T] extends GcsBatchIoCommand[T, StorageObject] {
  def file: GcsPath
  private val blob = file.blob
  override def operation: StorageRequest[StorageObject] = file.apiStorage.objects().get(blob.getBucket, blob.getName)
}

case class GcsBatchSizeCommand(override val file: GcsPath) extends IoSizeCommand(file) with GcsBatchGetCommand[Long] {
  override def mapGoogleResponse(response: StorageObject): Long = response.getSize.longValue()
}

case class GcsBatchCrc32Command(override val file: GcsPath) extends IoHashCommand(file) with GcsBatchGetCommand[String] {
  override def mapGoogleResponse(response: StorageObject): String = response.getCrc32c
}

case class GcsBatchTouchCommand(override val file: GcsPath) extends IoTouchCommand(file) with GcsBatchGetCommand[Unit] {
  override def mapGoogleResponse(response: StorageObject): Unit = ()
}