/*
 * Copyright 2018 Amazon.com, Inc. or its affiliates.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 *  BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 *  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 *  THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 *  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 *  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 *  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */
package cromwell.filesystems.s3.batch

import cromwell.core.callcaching.HashType.HashType
import cromwell.core.callcaching.{FileHashStrategy, HashType}
import software.amazon.awssdk.core.exception.SdkException
import software.amazon.awssdk.services.s3.model.{CopyObjectResponse, HeadObjectResponse, NoSuchKeyException}
import cromwell.core.io.{
  IoCommand,
  IoCopyCommand,
  IoDeleteCommand,
  IoExistsCommand,
  IoHashCommand,
  IoSizeCommand,
  IoTouchCommand
}
import cromwell.filesystems.s3.S3Path

/**
  * Io commands with S3 paths and some logic enabling batching of request.
  * @tparam T Return type of the IoCommand
  * @tparam U Return type of the response
  */
sealed trait S3BatchIoCommand[T, U] extends IoCommand[T] {

  /**
    * Maps the response of type U to the Cromwell Io response of type T
    */
  protected def mapResponse(response: U): T

  /**
    * Method called in the success callback of a batched request to decide what to do next.
    * Returns an `Either[T, S3BatchIoCommand[T, U]]`
    *   Left(value) means the command is complete, and the result can be sent back to the sender.
    *   Right(newCommand) means the command is not complete and needs another request to be executed.
    * Most commands will reply with Left(value).
    */
  def onSuccess(response: U): Either[T, S3BatchIoCommand[T, U]] =
    Left(mapResponse(response))

  /**
    * Override to handle a failure differently and potentially return a successful response.
    */
  def onFailure(error: SdkException): Option[Either[T, S3BatchIoCommand[T, U]]] = None
}

case class S3BatchCopyCommand(
  override val source: S3Path,
  override val destination: S3Path
) extends IoCopyCommand(source, destination)
    with S3BatchIoCommand[Unit, CopyObjectResponse] {
  override def mapResponse(response: CopyObjectResponse): Unit = ()
  override def commandDescription: String = s"S3BatchCopyCommand source '$source' destination '$destination'"
}

case class S3BatchDeleteCommand(
  override val file: S3Path,
  override val swallowIOExceptions: Boolean
) extends IoDeleteCommand(file, swallowIOExceptions)
    with S3BatchIoCommand[Unit, Void] {
  override protected def mapResponse(response: Void): Unit = ()
  override def commandDescription: String =
    s"S3BatchDeleteCommand file '$file' swallowIOExceptions '$swallowIOExceptions'"
}

/**
  * Base trait for commands that use the headObject() operation. (e.g: size, crc32, ...)
  */
sealed trait S3BatchHeadCommand[T] extends S3BatchIoCommand[T, HeadObjectResponse] {
  def file: S3Path
}

/**
  * `IoCommand` to find the size (content length) of an S3 object
  * @param file the path to the object
  */
case class S3BatchSizeCommand(override val file: S3Path) extends IoSizeCommand(file) with S3BatchHeadCommand[Long] {
  override def mapResponse(response: HeadObjectResponse): Long = response.contentLength
  override def commandDescription: String = s"S3BatchSizeCommand file '$file'"
}

/**
  * `IoCommand` to find the hash of an s3 object (the `Etag`)
  * @param file the path to the object
  */
case class S3BatchHashCommand(override val file: S3Path, override val hashStrategy: FileHashStrategy)
    extends IoHashCommand(file, hashStrategy)
    with S3BatchHeadCommand[String] {

  override def mapResponse(response: HeadObjectResponse): String =
    hashStrategy
      .getFileHash(
        response,
        (resp: HeadObjectResponse, hashType: HashType) =>
          hashType match {
            case HashType.Etag => Option(response.eTag())
            case _ => None
          }
      )
      .map(_.hash)
      .get
  // This ^^ unprotected .get is here because of a need to join the theoretical optionality of hashes
  // (codified in AN-380) with this class's lack of ability to handle optionality. Refactoring this class
  // isn't in scope right now. The .get replaced an unprotected response.eTag(), so we aren't creating any
  // new error cases. --jdewar 02-2025

  override def commandDescription: String = s"S3BatchEtagCommand file '$file' with hashStrategy '$hashStrategy'"
}

/**
  * `IoCommand` to "touch" an S3 object. The current implementation of `mapResponse` in this object doesn't do anything
  * as it is not clear that touch is meaningful in the context of S3
  * @param file the path to the object
  */
case class S3BatchTouchCommand(override val file: S3Path) extends IoTouchCommand(file) with S3BatchHeadCommand[Unit] {
  override def mapResponse(response: HeadObjectResponse): Unit = ()
  override def commandDescription: String = s"S3BatchTouchCommand file '$file'"
}

/**
  * `IoCommand` to determine the existence of an object in S3
  * @param file the path to the object
  */
case class S3BatchExistsCommand(override val file: S3Path)
    extends IoExistsCommand(file)
    with S3BatchHeadCommand[Boolean] {
  override def mapResponse(response: HeadObjectResponse): Boolean = true
  override def onFailure(error: SdkException): Option[Left[Boolean, Nothing]] =
    // If the object can't be found, don't fail the request but just return false as we were testing for existence
    error match {
      case _: NoSuchKeyException => Option(Left(false))
      case _ => None
    }
  override def commandDescription: String = s"S3BatchExistsCommand file '$file'"
}
