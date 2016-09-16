package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.{ActorRef, Props}
import better.files.File
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, FailedNonRetryableResponse, SucceededResponse}
import cromwell.backend.{BackendCacheHitCopyingActor, BackendJobDescriptor, BackendJobDescriptorKey}
import cromwell.core.logging.JobLogging
import cromwell.core.simpleton.{WdlValueBuilder, WdlValueSimpleton}
import cromwell.core.{JobOutputs, PathCopier}
import cromwell.services.metadata.CallMetadataKeys._
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata._
import wdl4s.values.WdlFile

import scala.util.{Failure, Success, Try}

case class JesCacheHitCopyingActor(override val jobDescriptor: BackendJobDescriptor,
                                   jesConfiguration: JesConfiguration,
                                   initializationData: JesBackendInitializationData,
                                   serviceRegistryActor: ActorRef) extends BackendCacheHitCopyingActor with JobLogging {

  override val configurationDescriptor = jesConfiguration.configurationDescriptor
  private lazy val runtimeAttributes = JesRuntimeAttributes(jobDescriptor.runtimeAttributes, jobLogger)
  private lazy val tag = s"${this.getClass.getSimpleName} [UUID(${workflowId.shortString}):${jobDescriptor.key.tag}]"
  private lazy val jesAttributes = jesConfiguration.jesAttributes
  private lazy val metadataJobKey = {
    val jobDescriptorKey: BackendJobDescriptorKey = jobDescriptor.key
    MetadataJobKey(jobDescriptorKey.call.fullyQualifiedName, jobDescriptorKey.index, jobDescriptorKey.attempt)
  }

  override def copyCachedOutputs(seqOfSimpletons: Seq[WdlValueSimpleton], jobDetritusFiles: Map[String,String],
                                 returnCode: Option[Int]): BackendJobExecutionResponse = {

    val gcsFileSystem = initializationData.workflowPaths.gcsFileSystemWithUserAuth
    val jesCallPaths = initializationData.workflowPaths.toJesCallPaths(jobDescriptor.key)
    val destinationCallRootPath: Path = jesCallPaths.callRootPath
    val getSourceCallRootPath: Try[Path] = Try {
      jobDetritusFiles.get(jesCallPaths.CallRootPathKey) match {
        case Some(srcPath) => gcsFileSystem.getPath(srcPath)
        case None => throw new RuntimeException(s"The call detritus files for source cache hit aren't found for call ${jobDescriptor.call.fullyQualifiedName}")
      }
    }

    def copyIfPossible: BackendJobExecutionResponse = {
      getSourceCallRootPath match {
        case Success(path: Path) =>
          copyAndRenameFiles(path)
          updateMetadata
          SucceededResponse(jobDescriptor.key, returnCode, createJobOutputs, Option(jobDetritusFiles))
        case Failure(error) => FailedNonRetryableResponse(jobDescriptor.key, error, None)
      }
    }

    def copyAndRenameFiles(sourceCallRootPath: Path) = {
      seqOfSimpletons foreach {
        case WdlValueSimpleton(key, wdlFile: WdlFile) =>
          copyFile(wdlFile.value, sourceCallRootPath)
        case _ => //ignoring non-WDLValueSimpleton:wdlFiles since they are not interesting here
      }

      jobDetritusFiles foreach {
        case (fileName: String, filePath: String) =>
          fileName match {
            case jesCallPaths.CallRootPathKey => // ignore
            case "stdout" => File(gcsFileSystem.getPath(filePath)).copyTo(jesCallPaths.stdoutPath, true)
            case "stderr" => File(gcsFileSystem.getPath(filePath)).copyTo(jesCallPaths.stderrPath, true)
            case "returnCode" => File(gcsFileSystem.getPath(filePath)).copyTo(jesCallPaths.returnCodePath, true)
            case "jesLog" => File(gcsFileSystem.getPath(filePath)).copyTo(jesCallPaths.jesLogPath, true)
            case "gcsExec" => File(gcsFileSystem.getPath(filePath)).copyTo(jesCallPaths.gcsExecPath, true)
          }
      }
      jobLogger.info("{} successfully finished copying all cached job outputs.", tag)
    }

    def copyFile(sourceString: String, srcContextPath: Path): Unit = {
      val sourcePath = gcsFileSystem.getPath(sourceString)
      PathCopier.copy(srcContextPath, sourcePath, destinationCallRootPath, overwrite = true)
    }

    def metadataEvent(key: String, value: Any) = {
      val metadataValue = MetadataValue(value)
      MetadataEvent(metadataKey(key), metadataValue)
    }

    def metadataKey(key: String) = MetadataKey(workflowId, Option(metadataJobKey), key)

    def updateMetadata = {
      val projectEvents = List(
        metadataEvent(JesMetadataKeys.GoogleProject, jesAttributes.project),
        metadataEvent(JesMetadataKeys.ExecutionBucket, jesAttributes.executionBucket),
        metadataEvent(JesMetadataKeys.EndpointUrl, jesAttributes.endpointUrl))

      val runtimeAttributesEvent = runtimeAttributes.asMap map {
        case (key, value) => MetadataEvent(metadataKey(s"runtimeAttributes:$key"), MetadataValue(value))
      }

      val detritusEvents = List(
          metadataEvent(CallMetadataKeys.CallRoot, destinationCallRootPath),
          metadataEvent(Stdout, jesCallPaths.stdoutPath),
          metadataEvent(Stderr, jesCallPaths.stderrPath),
          metadataEvent(ReturnCode, jesCallPaths.returnCodePath),
          metadataEvent(BackendLogsPrefix + ":log", jesCallPaths.jesLogPath),
          metadataEvent("cache:allowResultReuse", true)
          //TODO: Wire in copying of the monitoring log file some happy day
      )

      val events = runtimeAttributesEvent ++ projectEvents ++ detritusEvents

      serviceRegistryActor ! PutMetadataAction(events)
    }

    def createJobOutputs: JobOutputs = WdlValueBuilder.toJobOutputs(jobDescriptor.call.task.outputs, seqOfSimpletons)

    copyIfPossible
  }
}

object JesCacheHitCopyingActor {

  def props(jobDescriptor: BackendJobDescriptor,
            jesConfiguration: JesConfiguration,
            initializationData: JesBackendInitializationData,
            serviceRegistryActor: ActorRef): Props = {
    Props(new JesCacheHitCopyingActor(jobDescriptor, jesConfiguration, initializationData, serviceRegistryActor))
  }
}
