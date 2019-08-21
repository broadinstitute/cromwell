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

import com.google.cloud.storage.contrib.nio.CloudStorageOptions
import cromwell.backend.BackendInitializationData
import cromwell.backend.impl.aws.AwsBatchBackendInitializationData
import cromwell.backend.impl.aws.AWSBatchStorageSystems
import cromwell.backend.io.JobPaths
import cromwell.backend.standard.callcaching.{StandardCacheHitCopyingActor, StandardCacheHitCopyingActorParams}
import cromwell.core.io.DefaultIoCommandBuilder
import cromwell.core.CallOutputs
import cromwell.core.io.{IoCommand, IoTouchCommand}
import cromwell.core.path.Path
import cromwell.core.simpleton.{WomValueBuilder, WomValueSimpleton}
import cromwell.filesystems.s3.batch.S3BatchCommandBuilder
import common.util.TryUtil
import wom.values.WomFile

import scala.language.postfixOps
import scala.util.Try

class AwsBatchBackendCacheHitCopyingActor(standardParams: StandardCacheHitCopyingActorParams) extends StandardCacheHitCopyingActor(standardParams) {
  private val batchAttributes = BackendInitializationData
    .as[AwsBatchBackendInitializationData](standardParams.backendInitializationDataOption)
    .configuration.batchAttributes

  override protected val commandBuilder = batchAttributes.fileSystem match {
    case AWSBatchStorageSystems.s3  => S3BatchCommandBuilder
    case  _ => DefaultIoCommandBuilder
  }
  private val cachingStrategy = batchAttributes.duplicationStrategy

  override def processSimpletons(womValueSimpletons: Seq[WomValueSimpleton], sourceCallRootPath: Path) = (batchAttributes.fileSystem, cachingStrategy)  match {
    case (AWSBatchStorageSystems.s3 , UseOriginalCachedOutputs) =>  {
          val touchCommands: Seq[Try[IoTouchCommand]] = womValueSimpletons collect {
            case WomValueSimpleton(_, wdlFile: WomFile) => getPath(wdlFile.value) map S3BatchCommandBuilder.touchCommand
          }

          TryUtil.sequence(touchCommands) map {
            WomValueBuilder.toJobOutputs(jobDescriptor.taskCall.outputPorts, womValueSimpletons) -> _.toSet
          }
    }
    case (_, _) => super.processSimpletons(womValueSimpletons, sourceCallRootPath)
  }

  override def processDetritus(sourceJobDetritusFiles: Map[String, String]) = (batchAttributes.fileSystem, cachingStrategy)  match {
      case (AWSBatchStorageSystems.s3, UseOriginalCachedOutputs) => {
        // apply getPath on each detritus string file
        val detritusAsPaths = detritusFileKeys(sourceJobDetritusFiles).toSeq map { key =>
          key -> getPath(sourceJobDetritusFiles(key))
        } toMap

        // Don't forget to re-add the CallRootPathKey that has been filtered out by detritusFileKeys
        TryUtil.sequenceMap(detritusAsPaths, "Failed to make paths out of job detritus") map { newDetritus =>
          (newDetritus + (JobPaths.CallRootPathKey -> destinationCallRootPath)) -> newDetritus.values.map(S3BatchCommandBuilder.touchCommand).toSet
        }
      }
      case (_, _) => super.processDetritus(sourceJobDetritusFiles)
   }

  override protected def additionalIoCommands(sourceCallRootPath: Path,
                                              originalSimpletons: Seq[WomValueSimpleton],
                                              newOutputs: CallOutputs,
                                              originalDetritus:  Map[String, String],
                                              newDetritus: Map[String, Path]): List[Set[IoCommand[_]]] = (batchAttributes.fileSystem, cachingStrategy)  match {
      case (AWSBatchStorageSystems.s3, UseOriginalCachedOutputs) =>
            val content =
              s"""
                 |This directory does not contain any output files because this job matched an identical job that was previously run, thus it was a cache-hit.
                 |Cromwell is configured to not copy outputs during call caching. To change this, edit the filesystems.aws.caching.duplication-strategy field in your backend configuration.
                 |The original outputs can be found at this location: ${sourceCallRootPath.pathAsString}
      """.stripMargin

            List(Set(S3BatchCommandBuilder.writeCommand(jobPaths.callExecutionRoot / "call_caching_placeholder.txt", content, Seq(CloudStorageOptions.withMimeType("text/plain")))))
       case (AWSBatchStorageSystems.s3, CopyCachedOutputs) => List.empty
       case (_, _) => super.additionalIoCommands(sourceCallRootPath,originalSimpletons, newOutputs, originalDetritus,newDetritus)
    }
}
