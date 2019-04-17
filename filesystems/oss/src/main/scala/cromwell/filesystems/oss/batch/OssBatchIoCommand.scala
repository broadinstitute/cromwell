package cromwell.filesystems.oss.batch

import com.aliyun.oss.OSSException
import com.aliyun.oss.model._

import com.google.api.client.http.HttpHeaders

import cromwell.core.io._
import cromwell.filesystems.oss._

/**
  * Io commands with OSS paths and some logic enabling batching of request.
  * @tparam T Return type of the IoCommand
  * @tparam U Return type of the OSS response
  */
sealed trait OssBatchIoCommand[T, U] extends IoCommand[T] {
  /**
    * StorageRequest operation to be executed by this command
    */
  def operation: Any

  /**
    * Maps the Oss response of type U to the Cromwell Io response of type T
    */
  protected def mapOssResponse(response: U): T

  /**
    * Method called in the success callback of a batched request to decide what to do next.
    * Returns an Either[T, OssBatchIoCommand[T, U]]
    *   Left(value) means the command is complete, and the result can be sent back to the sender.
    *   Right(newCommand) means the command is not complete and needs another request to be executed.
    * Most commands will reply with Left(value).
    */
  def onSuccess(response: U, httpHeaders: HttpHeaders): Either[T, OssBatchIoCommand[T, U]] = {
    Left(mapOssResponse(response))
  }

  /**
    * Override to handle a failure differently and potentially return a successful response.
    */
  def onFailure(ossError: OSSException): Option[Either[T, OssBatchIoCommand[T, U]]] = None
}

case class OssBatchCopyCommand(
                                override val source: OssPath,
                                override val destination: OssPath,
                                override val overwrite: Boolean
                              ) extends IoCopyCommand(source, destination, overwrite) with OssBatchIoCommand[Unit, CopyObjectResult] {
  override def operation: GenericResult = {
    val getObjectRequest = new CopyObjectRequest(source.bucket, source.key, destination.bucket, destination.key)
    // TODO: Copy other attributes (encryption, metadata, etc.)
    source.ossClient.copyObject(getObjectRequest)
  }


  override def mapOssResponse(response: CopyObjectResult): Unit = ()
}

case class OssBatchDeleteCommand(
                                 override val file: OssPath,
                                 override val swallowIOExceptions: Boolean
                               ) extends IoDeleteCommand(file, swallowIOExceptions) with OssBatchIoCommand[Unit, Void] {
  def operation = file.ossClient.deleteObject(file.bucket, file.key)
  override protected def mapOssResponse(response: Void): Unit = ()
}

/**
  * Base trait for commands that use the headObject() operation. (e.g: size, crc32, ...)
  */
sealed trait OssBatchHeadCommand[T] extends OssBatchIoCommand[T, ObjectMetadata] {
  def file: OssPath

  override def operation: ObjectMetadata = file.ossClient.getObjectMetadata(file.bucket, file.key)
}

case class OssBatchSizeCommand(override val file: OssPath) extends IoSizeCommand(file) with OssBatchHeadCommand[Long] {
  override def mapOssResponse(response: ObjectMetadata): Long = response.getContentLength
}

case class OssBatchEtagCommand(override val file: OssPath) extends IoHashCommand(file) with OssBatchHeadCommand[String] {
  override def mapOssResponse(response: ObjectMetadata): String = response.getETag
}

case class OssBatchTouchCommand(override val file: OssPath) extends IoTouchCommand(file) with OssBatchHeadCommand[Unit] {
  override def mapOssResponse(response: ObjectMetadata): Unit = ()
}

case class OssBatchExistsCommand(override val file: OssPath) extends IoExistsCommand(file) with OssBatchIoCommand[Boolean, Boolean] {
  override def operation: Boolean = {
    file.ossClient.doesObjectExist(file.bucket, file.key)
  }

  override def mapOssResponse(response: Boolean): Boolean = response
}
