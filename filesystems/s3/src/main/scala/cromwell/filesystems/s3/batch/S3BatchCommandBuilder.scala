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

import cromwell.core.callcaching.FileHashStrategy
import cromwell.core.io.IoCommand.IOMetricsCallback
import cromwell.core.io.{IoCommandBuilder, PartialIoCommandBuilder}
import cromwell.core.path.Path
import cromwell.filesystems.s3.S3Path

import scala.util.Try

/**
  * Generates commands for IO operations on S3
  */
private case object PartialS3BatchCommandBuilder extends PartialIoCommandBuilder {
  override def sizeCommand: PartialFunction[Path, Try[S3BatchSizeCommand]] = { case path: S3Path =>
    Try(S3BatchSizeCommand(path))
  }

  override def deleteCommand: PartialFunction[(Path, Boolean), Try[S3BatchDeleteCommand]] = {
    case (path: S3Path, swallowIoExceptions) => Try(S3BatchDeleteCommand(path, swallowIoExceptions))
  }

  override def copyCommand: PartialFunction[(Path, Path), Try[S3BatchCopyCommand]] = {
    case (src: S3Path, dest: S3Path) => Try(S3BatchCopyCommand(src, dest))
  }

  override def hashCommand: PartialFunction[(Path, FileHashStrategy, IOMetricsCallback), Try[S3BatchHashCommand]] = {
    case (path: S3Path, s, _) =>
      Try(S3BatchHashCommand(path, s))
  }

  override def touchCommand: PartialFunction[Path, Try[S3BatchTouchCommand]] = { case path: S3Path =>
    Try(S3BatchTouchCommand(path))
  }

  override def existsCommand: PartialFunction[Path, Try[S3BatchExistsCommand]] = { case path: S3Path =>
    Try(S3BatchExistsCommand(path))
  }
}

/**
  * Builds `IoCommand` instances for specific actions on S3
  */
case object S3BatchCommandBuilder extends IoCommandBuilder(List(PartialS3BatchCommandBuilder))
