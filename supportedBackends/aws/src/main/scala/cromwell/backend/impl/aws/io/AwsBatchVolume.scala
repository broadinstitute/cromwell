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
import software.amazon.awssdk.services.batch.model.{Host, MountPoint, Volume}
import cromwell.core.path.{DefaultPathBuilder, Path}
import common.exception.MessageAggregation
import common.validation.ErrorOr._
import cromwell.backend.DiskPatterns
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
  // In AWS, disks are auto-sized so these patterns match simply "local-disk" or "/some/mnt"
  val MountedDiskPattern: Regex = raw"""^\s*(${DiskPatterns.Directory})\s*$$""".r
  val LocalDiskPattern: Regex = raw"""^\s*local-disk\s*$$""".r

  def parse(s: String): Try[AwsBatchVolume] = {

    val validation: ErrorOr[AwsBatchVolume] = s match {
      case LocalDiskPattern() =>
        Valid(AwsBatchWorkingDisk())
      case MountedDiskPattern(mountPoint) =>
        Valid(AwsBatchEmptyMountedDisk(DefaultPathBuilder.get(mountPoint)))
      // In addition to the AWS-specific patterns above, we can also fall back to PAPI-style patterns and ignore the size
      case DiskPatterns.WorkingDiskPattern(_, _) =>
        Valid(AwsBatchWorkingDisk())
      case DiskPatterns.MountedDiskPattern(mountPoint, _, fsType) =>
        Valid(AwsBatchEmptyMountedDisk(DefaultPathBuilder.get(mountPoint),fsType))
      case _ =>
        s"Disk strings should be of the format 'local-disk' or '/mount/point' but got: '$s'".invalidNel
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
  def fsType: String
  def getHostPath(id: Option[String]) : String =  {
    id match {
      case Some(id) => mountPoint.toAbsolutePath.pathAsString + "/" + id
      case None   =>   mountPoint.toAbsolutePath.pathAsString
    }
  }
  def toVolume(id: Option[String]=None): Volume = {
    Volume
      .builder
      .name(name)
      .host(Host.builder.sourcePath(getHostPath(id)).build)
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

case class AwsBatchEmptyMountedDisk(mountPoint: Path, ftype:String="ebs") extends AwsBatchVolume {
  val name = s"d-${mountPoint.pathAsString.md5Sum}"
  val fsType = ftype.toLowerCase
  override def toString: String = s"$name $mountPoint"
}

object AwsBatchWorkingDisk {
  val MountPoint: Path = DefaultPathBuilder.get("/cromwell_root")
  val Name = "local-disk"
  val fsType=  "ebs"
  val Default = AwsBatchWorkingDisk()
}

case class AwsBatchWorkingDisk() extends AwsBatchVolume {
  val mountPoint = AwsBatchWorkingDisk.MountPoint
  val name = AwsBatchWorkingDisk.Name
  val fsType = AwsBatchWorkingDisk.fsType
  override def toString: String = s"$name $mountPoint"
}
