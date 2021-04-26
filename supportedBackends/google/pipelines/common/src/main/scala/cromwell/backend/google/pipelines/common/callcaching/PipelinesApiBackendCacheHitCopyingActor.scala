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
import cromwell.core.WorkflowOptions.{CallCacheEgress, WorkflowOption}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import wom.values.WomFile

import scala.language.postfixOps
import scala.util.Try

class PipelinesApiBackendCacheHitCopyingActor(standardParams: StandardCacheHitCopyingActorParams) extends StandardCacheHitCopyingActor(standardParams) {
  override protected val commandBuilder: GcsBatchCommandBuilder.type = GcsBatchCommandBuilder
  private val cachingStrategy = BackendInitializationData
    .as[PipelinesApiBackendInitializationData](standardParams.backendInitializationDataOption)
    .papiConfiguration.papiAttributes.cacheHitDuplicationStrategy
  
  override def processSimpletons(womValueSimpletons: Seq[WomValueSimpleton],
                                 sourceCallRootPath: Path,
                                ): Try[(CallOutputs, Set[IoCommand[_]])] =
    cachingStrategy match {
    case CopyCachedOutputs =>
      // This is the super's implementation copy and pasted here, with an additional location check
      // If the location of the buckets cause an undesired egress charge, raise an exception
      val (destinationSimpletons, ioCommands): (List[WomValueSimpleton], Set[IoCommand[_]]) = womValueSimpletons.toList.foldMap({
        case WomValueSimpleton(key, wdlFile: WomSingleFile) =>
          val sourcePath = getPath(wdlFile.value).get
          val destinationPath = PathCopier.getDestinationFilePath(sourceCallRootPath, sourcePath, destinationCallRootPath)

          val sourceBucket = sourcePath.toString.split("/")(2)
          val destinationBucket = destinationPath.toString.split("/")(2)
          log.info(s"Comparing $sourceBucket and $destinationBucket")
          // IMPORTANT TODO: Obviously I can't just call bucket.location
          // How do I properly use GcsBatchCommandBuilder.locationCommand here?
          // And how do I get the WorkflowOption CallCacheEgress here?
          if (sourceBucket.location != destinationBucket.location) (failAndStop(CopyAttemptError(new TimeoutException("Buckets are from different locations"))))
          
          val destinationSimpleton = WomValueSimpleton(key, WomSingleFile(destinationPath.pathAsString))

          log.info(s"WILLY: copying $sourcePath to $destinationPath in processSimpletons")
          List(destinationSimpleton) -> Set(commandBuilder.copyCommand(sourcePath, destinationPath).get)
        case nonFileSimpleton => (List(nonFileSimpleton), Set.empty[IoCommand[_]])
      })
      log.info(s"WILLY: destinationSimpletons are $destinationSimpletons")
      log.info(s"WILLY: ioCommands are $ioCommands")
      (WomValueBuilder.toJobOutputs(jobDescriptor.taskCall.outputPorts, destinationSimpletons), ioCommands)
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
      val fileKeys = detritusFileKeys(sourceJobDetritusFiles)

      val zero = (Map.empty[String, Path], Set.empty[IoCommand[_]])

      val (destinationDetritus, ioCommands) = fileKeys.foldLeft(zero)({
        case ((detrituses, commands), detritus) =>
          val sourcePath = getPath(sourceJobDetritusFiles(detritus)).get
          val destinationPath = destinationJobDetritusPaths(detritus)

          val newDetrituses = detrituses + (detritus -> destinationPath)
        val sourceBucket = sourcePath.toString.split("/")(2)
        val destinationBucket = destinationPath.toString.split("/")(2)
        log.info(s"Comparing $sourceBucket and $destinationBucket")
        // IMPORTANT TODO: Obviously I can't just call bucket.location
        // How do I properly use GcsBatchCommandBuilder.locationCommand here?
        if (sourceBucket.location != destinationBucket.location) (failAndStop(CopyAttemptError(new TimeoutException("Buckets are from different locations"))))
        log.info(s"WILLY: copying $sourcePath to $destinationPath in processDetritus")
        (newDetrituses, commands + commandBuilder.copyCommand(sourcePath, destinationPath).get)
      })
      log.info(s"WILLY: destinationDetritus is $destinationDetritus")
      log.info(s"WILLY: ioCommands are $ioCommands")
      (destinationDetritus + (JobPaths.CallRootPathKey -> destinationCallRootPath), ioCommands)
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
