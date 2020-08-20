package cromwell.backend.google.pipelines.common.monitoring

import cromwell.backend.BackendJobDescriptor
import cromwell.backend.google.pipelines.common.WorkflowOptionKeys
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.MountsToEnv
import cromwell.backend.google.pipelines.common.io.PipelinesApiAttachedDisk
import cromwell.backend.io.WorkflowPaths
import cromwell.core.WorkflowOptions
import cromwell.core.path.{Path, PathFactory}

final class MonitoringImage(jobDescriptor: BackendJobDescriptor,
                            workflowOptions: WorkflowOptions,
                            workflowPaths: WorkflowPaths,
                            commandDirectory: Path,
                            workingDisk: PipelinesApiAttachedDisk,
                            localMonitoringImageScriptPath: Path,
                           ) {

  val monitoringImageOption: Option[String] = workflowOptions.get(WorkflowOptionKeys.MonitoringImage).toOption

  val monitoringImageScriptContainerPath: Path = workingDisk.mountPoint.resolve(localMonitoringImageScriptPath)

  val monitoringImageScriptOption: Option[Path] =
    for {
      _ <- monitoringImageOption // Only use the monitoring_image_script when monitoring_image provided
      monitoringImageScript <- workflowOptions.get(WorkflowOptionKeys.MonitoringImageScript).toOption
    } yield {
      PathFactory.buildPath(
        monitoringImageScript,
        workflowPaths.pathBuilders,
      )
    }

  val monitoringImageCommand: List[String] =
    monitoringImageScriptOption match {
      case Some(_) => List(
        "/bin/sh",
        "-c",
        s"cd '${commandDirectory.pathAsString}' && " +
          s"chmod +x '${monitoringImageScriptContainerPath.pathAsString}' && " +
          s"'${monitoringImageScriptContainerPath.pathAsString}'"
      )
      case None => Nil
    }

  val monitoringImageEnvironment: MountsToEnv = Env.monitoringImageEnvironment(jobDescriptor)
}
