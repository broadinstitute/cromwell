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
  override protected val commandBuilder: GcsBatchCommandBuilder.type = GcsBatchCommandBuilder
  private val cachingStrategy = BackendInitializationData
    .as[PipelinesApiBackendInitializationData](standardParams.backendInitializationDataOption)
    .papiConfiguration.papiAttributes.cacheHitDuplicationStrategy

  private def getBucketContinent(bucketLocation: String): String = {
    val indexOf = bucketLocation.indexOf("-")
    if (indexOf != -1) {
      bucketLocation.substring(0, indexOf)
    } else {
      bucketLocation
    }
  }

  override protected def locationCheckRequired: Boolean = {
    val option = standardParams.jobDescriptor.workflowDescriptor.workflowOptions.getOrElse("call_cache_egress", "global")
    val b = option != "global"
    log.warning(s"WILLY, in locationCheckRequired, returning $b because call_cache_egress option is $option")
    b
  }

  override protected def copyAllowedFromLocation(sourceLocation: String, destinationLocation: String): Boolean = {
    val callCacheEgressOption = standardParams.jobDescriptor.workflowDescriptor.workflowOptions.getOrElse("call_cache_egress", "global")
    log.warning(s"WILLY, in copyAllowedFromLocation, sourceLocation is $sourceLocation")
    log.warning(s"WILLY, in copyAllowedFromLocation, destinationLocation is $destinationLocation")
    log.warning(s"WILLY, in copyAllowedFromLocation, callCacheEgressOption is $callCacheEgressOption")
    callCacheEgressOption match {
      case "global" => true
      case "continental" =>
        val sourceContinent = getBucketContinent(sourceLocation)
        val destinationContinent = getBucketContinent(destinationLocation)
        val b = sourceContinent == destinationContinent
        log.warning(s"WILLY, in copyAllowedFromLocation, case continental, returning $b")
        b
      case "none" =>
        val b = sourceLocation == destinationLocation
        log.warning(s"WILLY, in copyAllowedFromLocation, case none, returning $b")
        b
      case default =>
        log.warning(s"WILLY, in copyAllowedFromLocation, did not match on callCacheEgressOption $default")
        true
    }    
  }

  override def processSimpletons(womValueSimpletons: Seq[WomValueSimpleton],
                                 sourceCallRootPath: Path,
                                ): Try[(CallOutputs, Set[IoCommand[_]])] =
    cachingStrategy match {
    case CopyCachedOutputs => super.processSimpletons(womValueSimpletons, sourceCallRootPath)
    case UseOriginalCachedOutputs =>
      val touchCommands: Seq[Try[IoTouchCommand]] = womValueSimpletons collect {
        case WomValueSimpleton(_, wdlFile: WomFile) => getPath(wdlFile.value) flatMap GcsBatchCommandBuilder.touchCommand
      }
      
      TryUtil.sequence(touchCommands) map {
        WomValueBuilder.toJobOutputs(jobDescriptor.taskCall.outputPorts, womValueSimpletons) -> _.toSet
      }
  }

  override def extractBlacklistPrefix(path: String): Option[String] = Option(path.stripPrefix("gs://").takeWhile(_ != '/'))

  override def processDetritus(sourceJobDetritusFiles: Map[String, String]
                              ): Try[(Map[String, Path], Set[IoCommand[_]])] = 
    cachingStrategy match {
    case CopyCachedOutputs => super.processDetritus(sourceJobDetritusFiles)
    case UseOriginalCachedOutputs =>
      // apply getPath on each detritus string file
      val detritusAsPaths = detritusFileKeys(sourceJobDetritusFiles).toSeq map { key =>
        key -> getPath(sourceJobDetritusFiles(key))
      } toMap

      // Don't forget to re-add the CallRootPathKey that has been filtered out by detritusFileKeys
      TryUtil.sequenceMap(detritusAsPaths, "Failed to make paths out of job detritus") flatMap { newDetritus =>
        Try {
          // PROD-444: Keep It Short and Simple: Throw on the first error and let the outer Try catch-and-re-wrap
          (newDetritus + (JobPaths.CallRootPathKey -> destinationCallRootPath)) ->
            newDetritus.values.map(GcsBatchCommandBuilder.touchCommand(_).get).toSet
        }
      }
  }

  override protected def additionalIoCommands(sourceCallRootPath: Path,
                                              originalSimpletons: Seq[WomValueSimpleton],
                                              newOutputs: CallOutputs,
                                              originalDetritus:  Map[String, String],
                                              newDetritus: Map[String, Path]): Try[List[Set[IoCommand[_]]]] = Try {
    cachingStrategy match {
      case UseOriginalCachedOutputs =>
        val content =
          s"""
             |This directory does not contain any output files because this job matched an identical job that was previously run, thus it was a cache-hit.
             |Cromwell is configured to not copy outputs during call caching. To change this, edit the filesystems.gcs.caching.duplication-strategy field in your backend configuration.
             |The original outputs can be found at this location: ${sourceCallRootPath.pathAsString}
      """.stripMargin

        // PROD-444: Keep It Short and Simple: Throw on the first error and let the outer Try catch-and-re-wrap
        List(Set(
          GcsBatchCommandBuilder.writeCommand(
            path = jobPaths.forCallCacheCopyAttempts.callExecutionRoot / "call_caching_placeholder.txt",
            content = content,
            options = Seq(CloudStorageOptions.withMimeType("text/plain")),
          ).get
        ))
      case CopyCachedOutputs => List.empty
    }
  }
}
