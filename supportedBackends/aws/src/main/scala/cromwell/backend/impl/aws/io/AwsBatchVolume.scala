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
package cromwell.backend.impl.aws.io

import cats.data.Validated._
import cats.syntax.validated._
import software.amazon.awssdk.services.batch.model.{MountPoint, Volume, Host}
import cromwell.core.path.{DefaultPathBuilder, Path}
import common.exception.MessageAggregation
import common.validation.ErrorOr._
import wom.values._

import scala.util.Try
import scala.util.matching.Regex


/*
 * This will handle volumes that are defined in the configuration. It will
 * *not* attach new block storage, as that is handled a priori as part of the
 * compute environment made available to run the jobs. This differs from
 * some other providers that create the entire compute environment on demand.
 */

object AwsBatchVolume {
  val Identifier = "[a-zA-Z0-9-_]{1,255}"
  val Directory = """/[^\s]+"""
  val WorkingDiskPattern: Regex = s"""${AwsBatchWorkingDisk.Name}\\s+($Identifier)""".r
  val MountedDiskPattern: Regex = s"""($Directory)\\s+($Identifier)""".r
  val LocalDiskPattern: Regex = s"""local-disk""".r

  def parse(s: String): Try[AwsBatchVolume] = {

    // TODO: This is probably significantly messed up and needs to be straightened
    // See BCS documentation here - we need volumes and mounts in a textual format
    // similar to what's defined: https://cromwell.readthedocs.io/en/develop/backends/BCS/
    // Thought: follow whatever format for the short form CLI version. Though
    // I'm not a fan of that format, it's at least consistent
    val validation: ErrorOr[AwsBatchVolume] = s match {
      case LocalDiskPattern() => Valid(AwsBatchWorkingDisk())
      case WorkingDiskPattern() => Valid(AwsBatchWorkingDisk())
      case MountedDiskPattern(mountPoint) => Valid(AwsBatchEmptyMountedDisk(DefaultPathBuilder.get(mountPoint)))
      case _ => s"Disk strings should be of the format 'local-disk' or '/mount/point SIZE TYPE' but got: '$s'".invalidNel
    }

    Try(validation match {
      case Valid(localDisk) => localDisk
      case Invalid(nels) =>
        throw new UnsupportedOperationException with MessageAggregation {
          val exceptionContext = ""
          val errorMessages: List[String] = nels.toList
        }
    })
  }
}

trait AwsBatchVolume {
  def name: String
  def mountPoint: Path
  def toVolume: Volume = {
    Volume
      .builder
      .name(name)
      .host(Host.builder.sourcePath(mountPoint.toAbsolutePath.pathAsString).build)
      .build
  }
  def toMountPoint: MountPoint = {
    MountPoint
      .builder
      .containerPath(mountPoint.toAbsolutePath.pathAsString)
      .sourceVolume(name)
      .build
  }
}

case class AwsBatchEmptyMountedDisk(mountPoint: Path) extends AwsBatchVolume {
  val name = s"d-${mountPoint.pathAsString.md5Sum}"
  override def toString: String = s"$name $mountPoint"
}

object AwsBatchWorkingDisk {
  val MountPoint: Path = DefaultPathBuilder.get("/cromwell_root")
  val Name = "local-disk"
  val Default = AwsBatchWorkingDisk()
}

case class AwsBatchWorkingDisk() extends AwsBatchVolume {
  val mountPoint = AwsBatchWorkingDisk.MountPoint
  val name = AwsBatchWorkingDisk.Name
  override def toString: String = s"$name $mountPoint"
}
