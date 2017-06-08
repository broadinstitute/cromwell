package cromwell.backend.impl.jes.callcaching

import cromwell.backend.BackendInitializationData
import cromwell.backend.impl.jes.JesBackendInitializationData
import cromwell.backend.standard.callcaching.{StandardCacheHitCopyingActor, StandardCacheHitCopyingActorParams}
import cromwell.core.path.Path
import cromwell.core.simpleton.{WdlValueBuilder, WdlValueSimpleton}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import lenthall.util.TryUtil

import scala.util.Success

class JesBackendCacheHitCopyingActor(standardParams: StandardCacheHitCopyingActorParams) extends StandardCacheHitCopyingActor(standardParams) with GcsBatchCommandBuilder {
  
  private val cachingStrategy = BackendInitializationData
    .as[JesBackendInitializationData](standardParams.backendInitializationDataOption)
    .jesConfiguration.jesAttributes.duplicationStrategy
  
  override def processSimpletons(wdlValueSimpletons: Seq[WdlValueSimpleton], sourceCallRootPath: Path) = cachingStrategy match {
    case CopyCachedOutputs => super.processSimpletons(wdlValueSimpletons, sourceCallRootPath)
    case UseOriginalCachedOutputs =>
      Success(WdlValueBuilder.toJobOutputs(jobDescriptor.call.task.outputs, wdlValueSimpletons) -> Set.empty[(Path, Path)])
  }
  
  override def processDetritus(sourceJobDetritusFiles: Map[String, String]) = cachingStrategy match {
    case CopyCachedOutputs => super.processDetritus(sourceJobDetritusFiles)
    case UseOriginalCachedOutputs =>
      TryUtil.sequenceMap(sourceJobDetritusFiles.mapValues(getPath), "Failed to make paths out of job detritus") map { _ -> Set.empty[(Path, Path)] }
  }
}
