package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, FailedNonRetryableResponse, SucceededResponse}
import cromwell.backend.{BackendCacheHitCopyingActor, BackendJobDescriptor, BackendJobDescriptorKey}
import cromwell.core.simpleton.{WdlValueBuilder, WdlValueSimpleton}
import cromwell.core.{JobOutputs, PathCopier}
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata._
import org.slf4j.LoggerFactory
import wdl4s.values.WdlFile

import scala.concurrent.Future

case class JesCacheHitCopyingActor(override val jobDescriptor: BackendJobDescriptor,
                                   jesConfiguration: JesConfiguration,
                                   initializationData: JesBackendInitializationData,
                                   serviceRegistryActor: ActorRef) extends BackendCacheHitCopyingActor {

  override val configurationDescriptor = jesConfiguration.configurationDescriptor
  private lazy val jesAttributes = jesConfiguration.jesAttributes
  private lazy val metadataJobKey = {
    val jobDescriptorKey: BackendJobDescriptorKey = jobDescriptor.key
    MetadataJobKey(jobDescriptorKey.call.fullyQualifiedName, jobDescriptorKey.index, jobDescriptorKey.attempt)
  }

  override def copyCachedOutputs(seqOfSimpletons: Seq[WdlValueSimpleton], jobDetritusFiles: Map[String,String],
                                 returnCode: Option[Int]): Future[BackendJobExecutionResponse] = {

      val gcsFileSystem = initializationData.workflowPaths.gcsFileSystemWithUserAuth
      val jesCallPaths = initializationData.workflowPaths.toJesCallPaths(jobDescriptor.key)
      val destinationCallRootPath: Path = jesCallPaths.callRootPath
      val getSourceCallRootPath = jobDetritusFiles.get(jesCallPaths.callRootPathKey) match {
        case Some(srcPath) => gcsFileSystem.getPath(srcPath)
        case None => val error = new RuntimeException(s"The call detritus files for source cache hit aren't found for call ${jobDescriptor.call.fullyQualifiedName}")
      }

    def cacheIfPossible = {
      getSourceCallRootPath match {
        case Some(path: Path) =>
          copyFiles(path)
          updateMetadata
          Future(SucceededResponse(jobDescriptor.key, returnCode, createJobOutputs))
        case error: RuntimeException => Future(FailedNonRetryableResponse(jobDescriptor.key, error, Option(-1)))
      }
    }

      def copyFiles(sourceCallRootPath: Path) = {
        seqOfSimpletons ++ jobDetritusFiles foreach {
          case WdlValueSimpleton(key, wdlFile: WdlFile) =>
            copyCachedCallOutputFiles(wdlFile.value, sourceCallRootPath)
          case (fileName: String, filePath: String) if fileName != jesCallPaths.callRootPathKey =>
            copyCachedCallOutputFiles(filePath, sourceCallRootPath)
          case _ =>
            val error = new RuntimeException(s"Unable to copy cached job outputs for ${jobDescriptor.call.fullyQualifiedName}")
            FailedNonRetryableResponse(jobDescriptor.key, error, Option(-1))
        }
      }

      def copyCachedCallOutputFiles(sourceString: String, srcContextPath: Path): Unit = {
        val sourcePath = gcsFileSystem.getPath(sourceString)
        PathCopier.copy(srcContextPath, sourcePath, destinationCallRootPath)
      }


      def tellMetadata(key: String, value: Any): Unit = {
        val event = metadataEvent(key, value)
        val putMetadataAction = PutMetadataAction(event)
        serviceRegistryActor ! putMetadataAction
      }

      def metadataEvent(key: String, value: Any) = {
        val metadataValue = MetadataValue(value)
        MetadataEvent(metadataKey(key), metadataValue)
      }

      def metadataKey(key: String) = MetadataKey(workflowId, Option(metadataJobKey), key)

      def updateMetadata = {
        tellMetadata(CallMetadataKeys.CallRoot, destinationCallRootPath)
        tellMetadata(JesMetadataKeys.GoogleProject, jesAttributes.project)
        tellMetadata(JesMetadataKeys.ExecutionBucket, jesAttributes.executionBucket)
        tellMetadata(JesMetadataKeys.EndpointUrl, jesAttributes.endpointUrl)
      }


    def createJobOutputs: JobOutputs = WdlValueBuilder.toJobOutputs(jobDescriptor.call.task.outputs, seqOfSimpletons)

    cacheIfPossible
  }
}

object JesCacheHitCopyingActor {
  val logger = LoggerFactory.getLogger("JesBackend")

  def props(jobDescriptor: BackendJobDescriptor,
            jesConfiguration: JesConfiguration,
            initializationData: JesBackendInitializationData,
            serviceRegistryActor: ActorRef): Props = {
    Props(new JesCacheHitCopyingActor(jobDescriptor, jesConfiguration, initializationData, serviceRegistryActor))
  }
}
