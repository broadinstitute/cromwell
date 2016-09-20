package cromwell.backend.callcaching

import java.nio.file.Path

import akka.actor.ActorRef
import cromwell.backend.BackendCacheHitCopyingActor
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, SucceededResponse}
import cromwell.backend.io.JobPaths
import cromwell.core.PathCopier
import cromwell.core.simpleton.{WdlValueBuilder, WdlValueSimpleton}
import wdl4s.values.WdlFile

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
  protected def getPath(file: String): Path

  protected def destinationCallRootPath: Path

  protected def serviceRegistryActor: ActorRef

  protected def destinationJobDetritusPaths: Map[String, Path]

  // Usually implemented by a subclass of JobCachingActorHelper
  protected def metadataKeyValues: Map[String, Any]

  private def lookupSourceCallRootPath(sourceJobDetritusFiles: Map[String, String]): Path = {
    sourceJobDetritusFiles.get(JobPaths.CallRootPathKey).map(getPath).getOrElse(throw new RuntimeException(
      s"The call detritus files for source cache hit aren't found for call ${jobDescriptor.call.fullyQualifiedName}")
    )
  }

  /**
    * After copying files, return the simpletons substituting the destination file paths.
    */
  private def copySimpletons(wdlValueSimpletons: Seq[WdlValueSimpleton],
                             sourceCallRootPath: Path): Seq[WdlValueSimpleton] = {
    wdlValueSimpletons map {
      case WdlValueSimpleton(key, wdlFile: WdlFile) =>
        val sourcePath = getPath(wdlFile.value)
        val destinationPath = PathCopier.getDestinationFilePath(sourceCallRootPath, sourcePath, destinationCallRootPath)
        duplicate(sourcePath, destinationPath)
        WdlValueSimpleton(key, WdlFile(destinationPath.toString))
      case wdlValueSimpleton => wdlValueSimpleton
    }
  }

  private def copyDetritus(sourceJobDetritusFiles: Map[String, String]): Map[String, String] = {
    val sourceKeys = sourceJobDetritusFiles.keySet
    val destinationKeys = destinationJobDetritusPaths.keySet
    val fileKeys = sourceKeys.intersect(destinationKeys).filterNot(_ == JobPaths.CallRootPathKey)

    val destinationJobDetritusFiles = fileKeys map { fileKey =>
      val sourcePath = getPath(sourceJobDetritusFiles(fileKey))
      val destinationPath = destinationJobDetritusPaths(fileKey)
      duplicate(sourcePath, destinationPath)
      (fileKey, destinationPath.toString)
    }

    destinationJobDetritusFiles.toMap
  }

  override def copyCachedOutputs(wdlValueSimpletons: Seq[WdlValueSimpleton],
                                 sourceJobDetritusFiles: Map[String, String],
                                 returnCodeOption: Option[Int]): BackendJobExecutionResponse = {
    val sourceCallRootPath = lookupSourceCallRootPath(sourceJobDetritusFiles)

    val destinationSimpletons = copySimpletons(wdlValueSimpletons, sourceCallRootPath)
    val destinationJobDetritusFiles = copyDetritus(sourceJobDetritusFiles)

    val destinationJobOutputs = WdlValueBuilder.toJobOutputs(jobDescriptor.call.task.outputs, destinationSimpletons)

    import cromwell.services.metadata.MetadataService.implicits.MetadataAutoPutter
    serviceRegistryActor.putMetadata(jobDescriptor.workflowDescriptor.id, Option(jobDescriptor.key), metadataKeyValues)

    SucceededResponse(jobDescriptor.key, returnCodeOption, destinationJobOutputs, Option(destinationJobDetritusFiles))
  }
}
