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

import common.util.TryUtil
import cromwell.backend.BackendInitializationData
import cromwell.backend.impl.aws.{
  AwsBatchBackendInitializationData,
  AwsBatchJobCachingActorHelper,
  AWSBatchStorageSystems
}
import cromwell.backend.io.JobPaths
// import cromwell.backend.impl.aws.io._

import cromwell.backend.standard.callcaching.{StandardCacheHitCopyingActor, StandardCacheHitCopyingActorParams}
import cromwell.core.CallOutputs
import cromwell.core.io.{DefaultIoCommandBuilder, IoCommand, IoCommandBuilder}
import cromwell.core.path.{Path, PathCopier} // ,DefaultPathBuilder}
import cromwell.core.simpleton.{WomValueBuilder, WomValueSimpleton}
import cromwell.filesystems.s3.batch.S3BatchCommandBuilder
import wom.values._ //WomFile
import wom.types._
import wom.callable.Callable.OutputDefinition
import wom.expression.NoIoFunctionSet
//import java.net.{URLDecoder}
//import cromwell.filesystems.s3.S3Path

import cats.syntax.validated._
import scala.language.postfixOps
import scala.util.Try //, Success}

class AwsBatchBackendCacheHitCopyingActor(standardParams: StandardCacheHitCopyingActorParams)
    extends StandardCacheHitCopyingActor(standardParams)
    with AwsBatchJobCachingActorHelper {
  private val batchAttributes = BackendInitializationData
    .as[AwsBatchBackendInitializationData](standardParams.backendInitializationDataOption)
    .configuration
    .batchAttributes

  override protected val commandBuilder: IoCommandBuilder = batchAttributes.fileSystem match {
    case AWSBatchStorageSystems.s3 => S3BatchCommandBuilder
    case _ => DefaultIoCommandBuilder
  }
  private val cachingStrategy = batchAttributes.duplicationStrategy

  // taken from the WomExpression trait:
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
  // taken from AwsBatchAsyncBackendJobExectuionHandler
  private def evaluateFiles(output: OutputDefinition): List[(WomFile, Boolean)] = {
    // mixed mandatory/optional types are cast to mandatory.
    val is_optional = areAllFilesOptional(output.womType)
    Try(
      output.expression
        .evaluateFiles(
          jobDescriptor.localInputs,
          NoIoFunctionSet,
          output.womType
        )
        .map(_.toList map { womFile =>
          // Pair each WomFile with the optional status
          (womFile.file, is_optional)
        })
    ).getOrElse(List.empty[(WomFile, Boolean)].validNel)
      .getOrElse(List.empty)
  }

  val womFileOutputsEvaluated = jobDescriptor.taskCall.callable.outputs
    .flatMap(evaluateFiles)

  // then get paths into a map of path => optional
  val womFileMap: Map[String, Boolean] = womFileOutputsEvaluated.map { case (womFile, isOptional) =>
    womFile.value -> isOptional
  }.toMap

  // starting from the womFileMap and simpletons, determine if the file is optional.
  // in womFileMap, the key is the path as speccified in the WDL (inside working dir), eg "new_dir/outfile.txt"
  // in WomValueSimpletons, the key is the full source path for cache copy, eg "s3://bucket/cromwell_temp/wf-id/call-id/new_dir/outfile.txt"
  // strategy :
  //    - check all keys in womFileMap. if simpleton ends with key, use optional_status from womFileMap (the value)
  //    - as soon as a mandatory file is found : break and return false
  //    - keep checking until the end otherwise, as multiple files can have nested suffixes (eg "new_dir/outfile.txt" and "outfile.txt")
  //    - if the key is found and optional : return true ; else return false
  def is_optional(womFile: String, womFileMap: Map[String, Boolean]): Boolean = {
    var isOptional = false
    for ((key, value) <- womFileMap)
      if (womFile.endsWith(key)) {
        if (!value) {
          return false
        }
        isOptional = value
      }
    return isOptional
  }

  // check if the file is on efs (local) or s3
  def is_efs(womFile: String): Boolean = {
    // get efs mount point/ disk from config
    val efs_mount = configuration.efsMntPoint.getOrElse("--")
    return womFile.startsWith(efs_mount)
  }

  override def processSimpletons(womValueSimpletons: Seq[WomValueSimpleton],
                                 sourceCallRootPath: Path
  ): Try[(CallOutputs, Set[IoCommand[_]])] =
    (batchAttributes.fileSystem, cachingStrategy) match {
      ///////////////////////
      // CACHE = REFERENCE //
      ///////////////////////
      case (AWSBatchStorageSystems.s3, UseOriginalCachedOutputs) =>
        val touchCommands: Seq[Try[(WomValueSimpleton, IoCommand[_])]] = womValueSimpletons collect {
          // only work on WomFiles
          case WomValueSimpleton(key, wdlFile: WomFile) =>
            val sourcePath = getPath(wdlFile.value).get
            // reference, so source == destination
            val destinationPath = sourcePath
            val destinationSimpleton = WomValueSimpleton(key, WomSingleFile(destinationPath.pathAsString))
            if (is_optional(wdlFile.value, womFileMap)) {
              // can I use this instead of noopCommand (from super) : case nonFileSimpleton => (List(nonFileSimpleton), Set.empty[IoCommand[_]])
              Try(destinationSimpleton -> S3BatchCommandBuilder.noopCommand(destinationPath).get)
            } else {
              Try(destinationSimpleton -> S3BatchCommandBuilder.existsOrThrowCommand(destinationPath).get)
            }
          case nonFileSimpleton =>
            Try(nonFileSimpleton -> S3BatchCommandBuilder.noopCommand(getPath("").get).get)
        }
        // group touchcommands
        TryUtil.sequence(touchCommands) map { simpletonsAndCommands =>
          val (destinationSimpletons, ioCommands) = simpletonsAndCommands.unzip
          WomValueBuilder.toJobOutputs(jobDescriptor.taskCall.outputPorts, destinationSimpletons) -> ioCommands.toSet
        }
      ///////////////////
      // CACHE == COPY //
      ///////////////////
      case (AWSBatchStorageSystems.s3, CopyCachedOutputs) =>
        val copyCommands: Seq[Try[(WomValueSimpleton, IoCommand[_])]] = womValueSimpletons collect {
          // only work on WomFiles
          case WomValueSimpleton(key, wdlFile: WomFile) =>
            val sourcePath = getPath(wdlFile.value).get
            // on efs : existOrthrow (mandatory) or noop (optional)
            if (is_efs(wdlFile.value)) {
              // on efs : source == destination
              val destinationPath = sourcePath
              val destinationSimpleton = WomValueSimpleton(key, WomSingleFile(destinationPath.pathAsString))
              if (is_optional(wdlFile.value, womFileMap)) {
                Try(destinationSimpleton -> S3BatchCommandBuilder.noopCommand(destinationPath).get)
              } else {
                Try(destinationSimpleton -> S3BatchCommandBuilder.existsOrThrowCommand(destinationPath).get)
              }
            }
            // on s3 : copy (mandatory) or copy if exists (optional)
            else {
              val destinationPath =
                PathCopier.getDestinationFilePath(sourceCallRootPath, sourcePath, destinationCallRootPath)
              val destinationSimpleton = WomValueSimpleton(key, WomSingleFile(destinationPath.pathAsString))

              // optional
              if (is_optional(wdlFile.value, womFileMap)) {
                val fileExists = sourcePath.exists
                if (fileExists) {
                  Try(destinationSimpleton -> S3BatchCommandBuilder.copyCommand(sourcePath, destinationPath).get)
                } else {
                  Try(destinationSimpleton -> S3BatchCommandBuilder.noopCommand(destinationPath).get)
                }

                // mandatory
              } else {
                Try(destinationSimpleton -> S3BatchCommandBuilder.copyCommand(sourcePath, destinationPath).get)
              }
            }
          case nonFileSimpleton =>
            Try(nonFileSimpleton -> S3BatchCommandBuilder.noopCommand(getPath("").get).get)
        }
        // get copycommands
        TryUtil.sequence(copyCommands) map { simpletonsAndCommands =>
          val (destinationSimpletons, ioCommands) = simpletonsAndCommands.unzip
          WomValueBuilder.toJobOutputs(jobDescriptor.taskCall.outputPorts, destinationSimpletons) -> ioCommands.toSet
        }
      ///////////////////////
      // NON-S3 FILESYSTEM //
      ///////////////////////
      case (_, _) => super.processSimpletons(womValueSimpletons, sourceCallRootPath)
    }

  // detritus files : job script, stdout, stderr and RC files.
  override def processDetritus(
    sourceJobDetritusFiles: Map[String, String]
  ): Try[(Map[String, Path], Set[IoCommand[_]])] =
    (batchAttributes.fileSystem, cachingStrategy) match {
      case (AWSBatchStorageSystems.s3, UseOriginalCachedOutputs) =>
        // apply getPath on each detritus string file
        val detritusAsPaths = detritusFileKeys(sourceJobDetritusFiles).toSeq map { key =>
          key -> getPath(sourceJobDetritusFiles(key))
        } toMap

        // Don't forget to re-add the CallRootPathKey that has been filtered out by detritusFileKeys
        TryUtil.sequenceMap(detritusAsPaths, "Failed to make paths out of job detritus") flatMap { newDetritus =>
          Try {
            // PROD-444: Keep It Short and Simple: Throw on the first error and let the outer Try catch-and-re-wrap
            (newDetritus + (JobPaths.CallRootPathKey -> destinationCallRootPath)) ->
              newDetritus.values.map(S3BatchCommandBuilder.existsOrThrowCommand(_).get).toSet
          }
        }
      case (_, _) => super.processDetritus(sourceJobDetritusFiles)
    }

  override protected def additionalIoCommands(sourceCallRootPath: Path,
                                              originalSimpletons: Seq[WomValueSimpleton],
                                              newOutputs: CallOutputs,
                                              originalDetritus: Map[String, String],
                                              newDetritus: Map[String, Path]
  ): Try[List[Set[IoCommand[_]]]] = Try {
    (batchAttributes.fileSystem, cachingStrategy) match {
      case (AWSBatchStorageSystems.s3, UseOriginalCachedOutputs) =>
        val content =
          s"""
             |This directory does not contain any output files because this job matched an identical job that was previously run, thus it was a cache-hit.
             |Cromwell is configured to not copy outputs during call caching. To change this, edit the filesystems.aws.caching.duplication-strategy field in your backend configuration.
             |The original outputs can be found at this location: ${sourceCallRootPath.pathAsString}
      """.stripMargin

        // PROD-444: Keep It Short and Simple: Throw on the first error and let the outer Try catch-and-re-wrap
        List(
          Set(
            S3BatchCommandBuilder
              .writeCommand(
                path = jobPaths.forCallCacheCopyAttempts.callExecutionRoot / "call_caching_placeholder.txt",
                content = content,
                options = Seq()
              )
              .get
          )
        )
      case (AWSBatchStorageSystems.s3, CopyCachedOutputs) => List.empty
      case (_, _) =>
        super
          .additionalIoCommands(sourceCallRootPath, originalSimpletons, newOutputs, originalDetritus, newDetritus)
          .get
    }
  }
}
