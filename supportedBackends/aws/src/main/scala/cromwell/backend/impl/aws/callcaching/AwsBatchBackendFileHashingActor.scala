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
package cromwell.backend.impl.aws.callcaching

import cromwell.backend.standard.callcaching.{StandardFileHashingActor, StandardFileHashingActorParams}
import cromwell.filesystems.s3.batch.S3BatchCommandBuilder
import cromwell.backend.BackendInitializationData
import cromwell.backend.impl.aws.AwsBatchBackendInitializationData
import cromwell.backend.impl.aws.AWSBatchStorageSystems
import cromwell.core.callcaching.FileHashStrategy
import cromwell.backend.impl.aws.AwsBatchJobCachingActorHelper
import cromwell.backend.io._
import cromwell.core.io.DefaultIoCommandBuilder
import common.validation.Validation._
import scala.util.Try
import cromwell.backend.standard.callcaching.StandardFileHashingActor.SingleFileHashRequest

import org.slf4j.{Logger, LoggerFactory}
import scala.util.Random
import wom.types.{
  WomArrayType,
  WomCompositeType,
  WomMapType,
  WomOptionalType,
  WomPairType,
  WomPrimitiveFileType,
  WomPrimitiveType,
  WomSingleFileType,
  WomType
}
import wom.values._

class AwsBatchBackendFileHashingActor(standardParams: StandardFileHashingActorParams)
    extends StandardFileHashingActor(standardParams)
    with AwsBatchJobCachingActorHelper {

  override val defaultHashingStrategies: Map[String, FileHashStrategy] = Map(
    ("s3", FileHashStrategy.ETag)
  )

  val Log: Logger = LoggerFactory.getLogger(StandardFileHashingActor.getClass)
  override val ioCommandBuilder = BackendInitializationData
    .as[AwsBatchBackendInitializationData](standardParams.backendInitializationDataOption)
    .configuration
    .batchAttributes
    .fileSystem match {
    case AWSBatchStorageSystems.s3 => S3BatchCommandBuilder
    case _ => DefaultIoCommandBuilder
  }
  // get backend config.
  val aws_config = BackendInitializationData
    .as[AwsBatchBackendInitializationData](standardParams.backendInitializationDataOption)
    .configuration

  // adapted from the WomExpression trait:
  private def areAllFilesOptional(womType: WomType): Boolean = {
    def innerAreAllFileTypesInWomTypeOptional(womType: WomType): Boolean = womType match {
      case WomOptionalType(_: WomPrimitiveFileType) =>
        true
      case _: WomPrimitiveFileType =>
        false
      case _: WomPrimitiveType =>
        true // WomPairTypes and WomCompositeTypes may have non-File components here which is fine.
      case WomArrayType(inner) => innerAreAllFileTypesInWomTypeOptional(inner)
      case WomMapType(_, inner) => innerAreAllFileTypesInWomTypeOptional(inner)
      case WomPairType(leftType, rightType) =>
        innerAreAllFileTypesInWomTypeOptional(leftType) && innerAreAllFileTypesInWomTypeOptional(rightType)
      case WomCompositeType(typeMap, _) => typeMap.values.forall(innerAreAllFileTypesInWomTypeOptional)
      case _ => false
    }

    // At the outermost level, primitives are never optional.
    womType match {
      case _: WomPrimitiveType =>
        false
      case _ => innerAreAllFileTypesInWomTypeOptional(womType)
    }
  }

  // custom strategy to handle optional and efs (local) files
  //  if optional file is missing : return hash of empty string (valid)
  //  if valid md5 is found for efs file: return the hash
  //  if no md5 is found : return None (pass request to parent hashing actor)
  //  if outdated md5 is found : return invalid string (assume file has been altered after md5 creation)
  //  if file is missing : return invalid string
  override def customHashStrategy(fileRequest: SingleFileHashRequest): Option[Try[String]] = {
    val file = getPath(fileRequest.file.value).get

    // the call inputs are mapped as path => optional
    val callInputFiles: Map[String, Boolean] = standardParams.jobDescriptor.fullyQualifiedInputs.flatMap {
      case (_, womFile) =>
        // Collect all WomArrays of WomFiles
        val arrays: Seq[WomArray] = womFile.collectAsSeq { case womFile: WomFile =>
          val files: List[WomSingleFile] = DirectoryFunctions
            .listWomSingleFiles(womFile, callPaths.workflowPaths)
            .toTry(s"Error getting single files for $womFile")
            .get
          WomArray(WomArrayType(WomSingleFileType), files)
        }
        // Determine if all files in the womFile are optional
        val isOptional = areAllFilesOptional(womFile.womType)
        // Flatten arrays and map each file path to its optional status
        arrays.flatMap(_.value).collect { case file: WomFile =>
          file.toString -> isOptional
        }
    }.toMap

    // optional files are allowed to be missing
    if (
      callInputFiles.contains(fileRequest.file.toString) && callInputFiles(fileRequest.file.toString) && !file.exists
    ) {
      // return hash of empty string
      Some("".md5Sum).map(str => Try(str))
      // the file is an efs file and sibling md5 is enabled
    } else if (
      file.toString.startsWith(aws_config.efsMntPoint.getOrElse("--")) && aws_config.checkSiblingMd5.getOrElse(false)
    ) {
      val md5 = file.sibling(s"${file.toString}.md5")
      // check existance of the main file :
      if (!file.exists) {
        // cache hit is invalid; return invalid (random) md5
        Some(Random.alphanumeric.take(32).mkString.md5Sum).map(str => Try(str))
      }
      // check existence of the sibling file and make sure it's newer than main file
      else if (md5.exists && md5.lastModifiedTime.isAfter(file.lastModifiedTime)) {
        // read the file.
        val md5_value: Option[String] = Some(md5.contentAsString.split("\\s+")(0))
        md5_value.map(str => Try(str))
      } else if (md5.exists && md5.lastModifiedTime.isBefore(file.lastModifiedTime)) {
        // sibling file is outdated, return invalid (random) string as md5
        Some(Random.alphanumeric.take(32).mkString.md5Sum).map(str => Try(str))
      } else {
        // File present, but no sibling found, fall back to default.
        None
      }
    } else {
      // non-efs file or sibling md5 is disabled : fall back to default
      None
    }
  }
}
