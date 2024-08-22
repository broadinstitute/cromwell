package cromwell.engine.workflow.lifecycle

import cromwell.backend.{AllBackendInitializationData, BackendConfigurationDescriptor, BackendInitializationData, BackendLifecycleActorFactory}
import cromwell.core.WorkflowOptions.UseRelativeOutputPaths
import cromwell.core.path.{Path, PathCopier, PathFactory}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.{BackendConfiguration, CromwellBackends}
import cromwell.engine.workflow.lifecycle.OutputsLocationHelper.FileRelocationMap
import wom.values.{WomSingleFile, WomValue}

object OutputsLocationHelper {
  type FileRelocationMap = Map[Path, Path]
}

trait OutputsLocationHelper {

  private def findFiles(values: Seq[WomValue]): Seq[WomSingleFile] =
    values flatMap {
      _.collectAsSeq { case file: WomSingleFile =>
        file
      }
    }

  protected def outputFilePathMapping(outputsDir: String,
                                      descriptor: EngineWorkflowDescriptor,
                                      backendInitData: AllBackendInitializationData,
                                      workflowOutputs: Seq[WomValue]
  ): FileRelocationMap = {
    val workflowOutputsPath = PathFactory.buildPath(outputsDir, descriptor.pathBuilders)
    val useRelativeOutputPaths: Boolean = descriptor.getWorkflowOption(UseRelativeOutputPaths).contains("true")
    val rootAndFiles = for {
      // NOTE: Without .toSeq, outputs in arrays only yield the last output
      backend <- descriptor.backendAssignments.values.toSeq
      config <- BackendConfiguration.backendConfigurationDescriptor(backend).toOption.toSeq
      rootPath <- getBackendRootPath(backend, config, descriptor, backendInitData).toSeq
      outputFiles = findFiles(workflowOutputs).map(_.value)
    } yield (rootPath, outputFiles)

    // This regex will make sure the path is relative to the execution folder.
    // the call-.* part is there to prevent arbitrary folders called execution to get caught.
    // Truncate regex is declared here. If it were declared in the if statement the regex would have to be
    // compiled for every single file.
    // "execution" should be optional, because its not created on AWS.
    // Also cacheCopy or attempt-<int> folders are optional.
    lazy val truncateRegex = ".*/call-[^/]*/(shard-[0-9]+/)?(cacheCopy/)?(attempt-[0-9]+/)?(execution/)?".r
    val outputFileDestinations = rootAndFiles flatMap { case (workflowRoot, outputs) =>
      outputs map { output =>
        val outputPath = PathFactory.buildPath(output, descriptor.pathBuilders)
        outputPath -> {
          if (useRelativeOutputPaths) {
            val pathRelativeToExecDir = truncateRegex.replaceFirstIn(outputPath.pathAsString, "")
            workflowOutputsPath.resolve(pathRelativeToExecDir)
          } else PathCopier.getDestinationFilePath(workflowRoot, outputPath, workflowOutputsPath)
        }
      }
    }
    outputFileDestinations.distinct.toMap
  }

  private def getBackendRootPath(backend: String,
                                 config: BackendConfigurationDescriptor,
                                 descriptor: EngineWorkflowDescriptor,
                                 backendInitData: AllBackendInitializationData
  ): Option[Path] =
    getBackendFactory(backend) map getRootPath(config, backendInitData.get(backend), descriptor)

  private def getBackendFactory(backend: String): Option[BackendLifecycleActorFactory] =
    CromwellBackends.backendLifecycleFactoryActorByName(backend).toOption

  private def getRootPath(config: BackendConfigurationDescriptor,
                          initializationData: Option[BackendInitializationData],
                          descriptor: EngineWorkflowDescriptor
  )(backendFactory: BackendLifecycleActorFactory): Path =
    backendFactory.getExecutionRootPath(descriptor.backendDescriptor, config.backendConfig, initializationData)

}
