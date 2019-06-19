package cromwell.backend.google.pipelines.common.callcaching

import com.google.cloud.storage.contrib.nio.CloudStorageOptions
import common.util.TryUtil
import cromwell.backend.BackendInitializationData
import cromwell.backend.google.pipelines.common.PipelinesApiBackendInitializationData
import cromwell.backend.io.JobPaths
import cromwell.backend.standard.callcaching.{StandardCacheHitCopyingActor, StandardCacheHitCopyingActorParams}
import cromwell.core.CallOutputs
import cromwell.core.io.{IoCommand, IoTouchCommand}
import cromwell.core.path.Path
import cromwell.core.simpleton.{WomValueBuilder, WomValueSimpleton}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import wom.values.WomFile

import scala.language.postfixOps
import scala.util.Try

class PipelinesApiBackendCacheHitCopyingActor(standardParams: StandardCacheHitCopyingActorParams) extends StandardCacheHitCopyingActor(standardParams) {
  override protected val commandBuilder = GcsBatchCommandBuilder
  private val cachingStrategy = BackendInitializationData
    .as[PipelinesApiBackendInitializationData](standardParams.backendInitializationDataOption)
    .papiConfiguration.papiAttributes.cacheHitDuplicationStrategy
  
  override def processSimpletons(womValueSimpletons: Seq[WomValueSimpleton], sourceCallRootPath: Path) = cachingStrategy match {
    case CopyCachedOutputs => super.processSimpletons(womValueSimpletons, sourceCallRootPath)
    case UseOriginalCachedOutputs =>
      val touchCommands: Seq[Try[IoTouchCommand]] = womValueSimpletons collect {
        case WomValueSimpleton(_, wdlFile: WomFile) => getPath(wdlFile.value) map GcsBatchCommandBuilder.touchCommand
      }
      
      TryUtil.sequence(touchCommands) map {
        WomValueBuilder.toJobOutputs(jobDescriptor.taskCall.outputPorts, womValueSimpletons) -> _.toSet
      }
  }

  override def extractBlacklistPrefix(path: String): Option[String] = Option(path.stripPrefix("gs://").takeWhile(_ != '/'))

  override def processDetritus(sourceJobDetritusFiles: Map[String, String]) = cachingStrategy match {
    case CopyCachedOutputs => super.processDetritus(sourceJobDetritusFiles)
    case UseOriginalCachedOutputs =>
      // apply getPath on each detritus string file
      val detritusAsPaths = detritusFileKeys(sourceJobDetritusFiles).toSeq map { key =>
        key -> getPath(sourceJobDetritusFiles(key))
      } toMap

      // Don't forget to re-add the CallRootPathKey that has been filtered out by detritusFileKeys
      TryUtil.sequenceMap(detritusAsPaths, "Failed to make paths out of job detritus") map { newDetritus =>
        (newDetritus + (JobPaths.CallRootPathKey -> destinationCallRootPath)) -> newDetritus.values.map(GcsBatchCommandBuilder.touchCommand).toSet
      }
  }

  override protected def additionalIoCommands(sourceCallRootPath: Path,
                                              originalSimpletons: Seq[WomValueSimpleton],
                                              newOutputs: CallOutputs,
                                              originalDetritus:  Map[String, String],
                                              newDetritus: Map[String, Path]): List[Set[IoCommand[_]]] = {
    cachingStrategy match {
      case UseOriginalCachedOutputs =>
        val content =
          s"""
             |This directory does not contain any output files because this job matched an identical job that was previously run, thus it was a cache-hit.
             |Cromwell is configured to not copy outputs during call caching. To change this, edit the filesystems.gcs.caching.duplication-strategy field in your backend configuration.
             |The original outputs can be found at this location: ${sourceCallRootPath.pathAsString}
      """.stripMargin

        List(Set(GcsBatchCommandBuilder.writeCommand(jobPaths.callExecutionRoot / "call_caching_placeholder.txt", content, Seq(CloudStorageOptions.withMimeType("text/plain")))))
      case CopyCachedOutputs => List.empty
    }
  }
}
