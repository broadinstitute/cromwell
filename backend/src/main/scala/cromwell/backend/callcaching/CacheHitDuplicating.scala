package cromwell.backend.callcaching

import java.nio.file.Path

import akka.actor.ActorRef
import cromwell.backend.BackendCacheHitCopyingActor
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, JobSucceededResponse}
import cromwell.backend.io.JobPaths
import cromwell.core.path.PathCopier
import cromwell.core.path.PathImplicits._
import cromwell.core.simpleton.{WdlValueBuilder, WdlValueSimpleton}
import wdl4s.values.WdlFile

import scala.language.postfixOps
import scala.util.Try

/**
  * Mixin implementing common functionality for a BackendCacheHitCopyingActor.
  *
  * Implements copyCachedOutputs, with abstract methods for converting a string to a path, duplicating a path, returning
  * a reference to the service registry actor, and obtaining various metadata/outputs for the job.
  */
trait CacheHitDuplicating {
  this: BackendCacheHitCopyingActor =>

  /**
    * Duplicates two paths returned by getPath.
    *
    * @param source      Source path.
    * @param destination Destination path.
    */
  protected def duplicate(source: Path, destination: Path): Unit

  /**
    * Returns an absolute path to the file.
    *
    * NOTE: If necessary for separate credentialing, we may split this method into getSourcePath and getDestinationPath.
    *
    * @param file the string version of the path
    * @return an absolute path to the file with potential credentials embedded within.
    */
  protected def getPath(file: String): Try[Path]

  protected def destinationCallRootPath: Path

  protected def serviceRegistryActor: ActorRef

  protected def destinationJobDetritusPaths: Map[String, Path]

  // Usually implemented by a subclass of JobCachingActorHelper
  protected def startMetadataKeyValues: Map[String, Any]

  private def lookupSourceCallRootPath(sourceJobDetritusFiles: Map[String, String]): Path = {
    sourceJobDetritusFiles.get(JobPaths.CallRootPathKey).map(getPath).get recover {
      case failure =>
        throw new RuntimeException(s"${JobPaths.CallRootPathKey} wasn't found for call ${jobDescriptor.call.fullyQualifiedName}", failure)
    } get
  }

  /**
    * After copying files, return the simpletons substituting the destination file paths.
    */
  private def copySimpletons(wdlValueSimpletons: Seq[WdlValueSimpleton],
                             sourceCallRootPath: Path): Seq[WdlValueSimpleton] = {
    wdlValueSimpletons map {
      case WdlValueSimpleton(key, wdlFile: WdlFile) =>
        val sourcePath = getPath(wdlFile.value).get
        val destinationPath = PathCopier.getDestinationFilePath(sourceCallRootPath, sourcePath, destinationCallRootPath)
        duplicate(sourcePath, destinationPath)
        WdlValueSimpleton(key, WdlFile(destinationPath.toRealString))
      case wdlValueSimpleton => wdlValueSimpleton
    }
  }

  private def copyDetritus(sourceJobDetritusFiles: Map[String, String]): Map[String, Path] = {
    val sourceKeys = sourceJobDetritusFiles.keySet
    val destinationKeys = destinationJobDetritusPaths.keySet
    val fileKeys = sourceKeys.intersect(destinationKeys).filterNot(_ == JobPaths.CallRootPathKey)

    val destinationJobDetritusFiles = fileKeys map { fileKey =>
      val sourcePath = getPath(sourceJobDetritusFiles(fileKey)).get
      val destinationPath = destinationJobDetritusPaths(fileKey)
      duplicate(sourcePath, destinationPath)
      (fileKey, destinationPath)
    }

    destinationJobDetritusFiles.toMap + (JobPaths.CallRootPathKey -> destinationCallRootPath)
  }

  override def copyCachedOutputs(wdlValueSimpletons: Seq[WdlValueSimpleton],
                                 sourceJobDetritusFiles: Map[String, String],
                                 returnCodeOption: Option[Int]): BackendJobExecutionResponse = {
    val sourceCallRootPath = lookupSourceCallRootPath(sourceJobDetritusFiles)

    val destinationSimpletons = copySimpletons(wdlValueSimpletons, sourceCallRootPath)
    val destinationJobDetritusFiles = copyDetritus(sourceJobDetritusFiles)

    val destinationJobOutputs = WdlValueBuilder.toJobOutputs(jobDescriptor.call.task.outputs, destinationSimpletons)

    import cromwell.services.metadata.MetadataService.implicits.MetadataAutoPutter
    serviceRegistryActor.putMetadata(
      jobDescriptor.workflowDescriptor.id, Option(jobDescriptor.key), startMetadataKeyValues)

    JobSucceededResponse(jobDescriptor.key, returnCodeOption, destinationJobOutputs, Option(destinationJobDetritusFiles), Seq.empty)
  }
}
