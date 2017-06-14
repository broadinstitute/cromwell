package cromwell.backend.impl.jes.callcaching

import com.google.cloud.storage.contrib.nio.CloudStorageOptions
import cromwell.backend.BackendInitializationData
import cromwell.backend.impl.jes.JesBackendInitializationData
import cromwell.backend.io.JobPaths
import cromwell.backend.standard.callcaching.{StandardCacheHitCopyingActor, StandardCacheHitCopyingActorParams}
import cromwell.core.CallOutputs
import cromwell.core.io.IoCommand
import cromwell.core.path.Path
import cromwell.core.simpleton.{WdlValueBuilder, WdlValueSimpleton}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import lenthall.util.TryUtil
import wdl4s.values.WdlFile

import scala.util.Success

class JesBackendCacheHitCopyingActor(standardParams: StandardCacheHitCopyingActorParams) extends StandardCacheHitCopyingActor(standardParams) with GcsBatchCommandBuilder {
  
  private val cachingStrategy = BackendInitializationData
    .as[JesBackendInitializationData](standardParams.backendInitializationDataOption)
    .jesConfiguration.jesAttributes.duplicationStrategy
  
  override def processSimpletons(wdlValueSimpletons: Seq[WdlValueSimpleton], sourceCallRootPath: Path) = cachingStrategy match {
    case CopyCachedOutputs => super.processSimpletons(wdlValueSimpletons, sourceCallRootPath)
    case UseOriginalCachedOutputs =>
      val touchCommands = wdlValueSimpletons collect {
        case WdlValueSimpleton(_, wdlFile: WdlFile) => touchCommand(getPath(wdlFile.value).get)
      }
      Success(WdlValueBuilder.toJobOutputs(jobDescriptor.call.task.outputs, wdlValueSimpletons) -> touchCommands.toSet)
  }
  
  override def processDetritus(sourceJobDetritusFiles: Map[String, String]) = cachingStrategy match {
    case CopyCachedOutputs => super.processDetritus(sourceJobDetritusFiles)
    case UseOriginalCachedOutputs =>
      val touchCommands: Set[IoCommand[_]] = detritusFileKeys(sourceJobDetritusFiles) map { key =>
        touchCommand(getPath(sourceJobDetritusFiles(key)).get)
      }
      TryUtil.sequenceMap(sourceJobDetritusFiles.mapValues(getPath), "Failed to make paths out of job detritus") map { _ -> touchCommands }
  }

  override protected def additionalIoCommands(newOutputs: CallOutputs, newDetritus: Map[String, Path]): List[Set[IoCommand[_]]] = {
    val originalExecutionRoot = newDetritus(JobPaths.CallRootPathKey)
    val content =
      s"""
         |This directory does not contain any output files because this job matched an identical job that was previously run, thus it was a cache-hit.
         |Cromwell is configured to not copy outputs during call caching. To change this, edit the filesystems.gcs.caching.duplication-strategy field in your backend configuration.
         |The original outputs can be found at this location: ${originalExecutionRoot.pathAsString}
      """.stripMargin

    List(Set(writeCommand(jobPaths.callExecutionRoot / "call_caching_placeholder.txt", content, Seq(CloudStorageOptions.withMimeType("text/plain")))))
  }
}
