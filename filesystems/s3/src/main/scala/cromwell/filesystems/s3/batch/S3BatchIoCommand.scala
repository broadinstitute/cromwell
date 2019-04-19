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

import software.amazon.awssdk.core.exception.SdkException

import software.amazon.awssdk.services.s3.model.{HeadObjectResponse, CopyObjectResponse, NoSuchKeyException}
import cromwell.core.io.{IoCommand,
                         IoDeleteCommand,
                         IoSizeCommand,
                         IoHashCommand,
                         IoTouchCommand,
                         IoExistsCommand,
                         IoCopyCommand}

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
    * Returns an Either[T, S3BatchIoCommand[T, U]]
    *   Left(value) means the command is complete, and the result can be sent back to the sender.
    *   Right(newCommand) means the command is not complete and needs another request to be executed.
    * Most commands will reply with Left(value).
    */
  def onSuccess(response: U): Either[T, S3BatchIoCommand[T, U]] = {
    Left(mapResponse(response))
  }

  /**
    * Override to handle a failure differently and potentially return a successful response.
    */
  def onFailure(error: SdkException): Option[Either[T, S3BatchIoCommand[T, U]]] = None
}

case class S3BatchCopyCommand(
                           override val source: S3Path,
                           override val destination: S3Path,
                           override val overwrite: Boolean,
                         ) extends IoCopyCommand(source, destination, overwrite) with S3BatchIoCommand[Unit, CopyObjectResponse] {
  override def mapResponse(response: CopyObjectResponse): Unit = ()
}

case class S3BatchDeleteCommand(
                                  override val file: S3Path,
                                  override val swallowIOExceptions: Boolean
                                ) extends IoDeleteCommand(file, swallowIOExceptions) with S3BatchIoCommand[Unit, Void] {
  override protected def mapResponse(response: Void): Unit = ()
}

/**
  * Base trait for commands that use the headObject() operation. (e.g: size, crc32, ...)
  */
sealed trait S3BatchHeadCommand[T] extends S3BatchIoCommand[T, HeadObjectResponse] {
  def file: S3Path
}

case class S3BatchSizeCommand(override val file: S3Path) extends IoSizeCommand(file) with S3BatchHeadCommand[Long] {
  override def mapResponse(response: HeadObjectResponse): Long = response.contentLength
}

case class S3BatchEtagCommand(override val file: S3Path) extends IoHashCommand(file) with S3BatchHeadCommand[String] {
  override def mapResponse(response: HeadObjectResponse): String = response.eTag
}

case class S3BatchTouchCommand(override val file: S3Path) extends IoTouchCommand(file) with S3BatchHeadCommand[Unit] {
  override def mapResponse(response: HeadObjectResponse): Unit = ()
}

case class S3BatchExistsCommand(override val file: S3Path) extends IoExistsCommand(file) with S3BatchHeadCommand[Boolean] {
  override def mapResponse(response: HeadObjectResponse): Boolean = true

  override def onFailure(error: SdkException) = {
    // If the object can't be found, don't fail the request but just return false as we were testing for existence
    error match {
      case _ : NoSuchKeyException => Option(Left(false))
      case _ => None
    }
  }
}
