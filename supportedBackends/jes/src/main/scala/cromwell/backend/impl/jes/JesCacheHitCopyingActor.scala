package cromwell.backend.impl.jes

import akka.actor.{ActorRef, Props}
import cromwell.backend.{BackendCacheHitCopyingActor, BackendJobDescriptor}
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, SucceededResponse}
import cromwell.backend.impl.jes.JesImplicits.PathString
import cromwell.core.{JobOutputs, PathCopier}
import org.slf4j.LoggerFactory

import scala.concurrent.Future

object JesCacheHitCopyingActor {
  val logger = LoggerFactory.getLogger("JesBackend")

  def props(jobDescriptor: BackendJobDescriptor,
            jesConfiguration: JesConfiguration,
            initializationData: JesBackendInitializationData,
            serviceRegistryActor: ActorRef): Props = {
    Props(new JesCacheHitCopyingActor(jobDescriptor, jesConfiguration, initializationData, serviceRegistryActor))
  }
}

case class JesCacheHitCopyingActor(override val jobDescriptor: BackendJobDescriptor,
                                   jesConfiguration: JesConfiguration,
                                   initializationData: JesBackendInitializationData,
                                   serviceRegistryActor: ActorRef) extends BackendCacheHitCopyingActor {

  override val configurationDescriptor = jesConfiguration.configurationDescriptor

  override def copyCachedOutputs(cachedJobOutputs: JobOutputs): Future[BackendJobExecutionResponse] = {
    val workflowPaths = initializationData.workflowPaths
    val jesWorkflowPaths = workflowPaths.toJesCallPaths(jobDescriptor.key)
    val callRootPath = jesWorkflowPaths.callRootPath
    val callRootParent = callRootPath.getParent
    val callRootRoot = callRootPath.getRoot

    //renaming requires some tech talk
    def renameCallSpecificFiles: List[String] = {
      val jesFileNames = List(
        jesWorkflowPaths.returnCodeFilename,
        jesWorkflowPaths.stderrFilename,
        jesWorkflowPaths.stdoutFilename,
        jesWorkflowPaths.jesLogFilename
      )
      jesFileNames
    }

    val gcsFileSystem = workflowPaths.gcsFileSystemWithUserAuth

    def copyFilesInGcs() = {
      val stringFilePaths = cachedJobOutputs map { case (lqn, jobOutput) => lqn -> jobOutput.wdlValue.valueString }

      stringFilePaths foreach {
        case (lqn, filePath) => if (filePath.isGcsUrl)
          copyCachedCallOutputFiles(filePath)
      }
    }

    //nested array output copying? or more complicated file structures
    def copyCachedCallOutputFiles(sourcePath: String): Unit = {
      val sourceCallDirectory = gcsFileSystem.getPath(sourcePath)
      // suspicious that this is unused
      val buildCallDestinationDirectory = gcsFileSystem.getPathAsDirectory(callRootPath.getParent.toString)

      PathCopier.copy(callRootPath, sourceCallDirectory.getParent, callRootParent)
    }

    copyFilesInGcs()
    Future(SucceededResponse(jobDescriptor.key, Option(0), cachedJobOutputs))
  }
}
