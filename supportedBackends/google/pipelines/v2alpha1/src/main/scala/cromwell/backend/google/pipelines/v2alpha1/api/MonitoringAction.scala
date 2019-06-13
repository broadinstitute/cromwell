package cromwell.backend.google.pipelines.v2alpha1.api

import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import cromwell.backend.google.pipelines.common.WorkflowOptionKeys
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.io.PipelinesApiAttachedDisk

import scala.collection.JavaConverters._
import spray.json._
import spray.json.DefaultJsonProtocol._
import wom.values.{WomFile, WomValue}

trait MonitoringAction {
  object Env {
    /**
      * Name of an environmental variable
      */
    val WorkflowId = "WORKFLOW_ID"
    val WorkflowName = "WORKFLOW_NAME"
    val TaskCallName = "TASK_CALL_NAME"
    val TaskCallIndex = "TASK_CALL_INDEX"
    val TaskCallAttempt = "TASK_CALL_ATTEMPT"
    val TaskInputs = "TASK_INPUTS"
    val TaskDisks = "TASK_DISKS"
    val DiskMounts = "DISK_MOUNTS" // preserved for backward compatibility
    val MonitoringConfig = "MONITORING_CONFIG"
  }

  private def monitoringAction(createPipelineParameters: CreatePipelineParameters, image: String, config: String, mounts: List[Mount]): List[Action] = {
    val job = createPipelineParameters.jobDescriptor

    val environment = Map(
      Env.WorkflowId -> job.workflowDescriptor.id.toString,
      Env.WorkflowName -> job.workflowDescriptor.callable.name,
      Env.TaskCallName -> job.taskCall.localName,
      Env.TaskCallIndex -> (job.key.index map { _.toString } getOrElse "NA"),
      Env.TaskCallAttempt -> job.key.attempt.toString,
      Env.TaskInputs -> getTaskInputs(job.fullyQualifiedInputs),
      Env.TaskDisks -> getTaskDisks(createPipelineParameters.adjustedSizeDisks, mounts),
      Env.DiskMounts -> mounts.map(_.getPath).mkString(" "),
      Env.MonitoringConfig -> config,
    )

    val monitoringAction = ActionBuilder.monitoringAction(image, environment, mounts)

    val describeAction = ActionBuilder.describeDocker("monitoring action", monitoringAction)
      .setFlags(List(ActionFlag.RunInBackground.toString).asJava)

    val terminationAction = ActionBuilder.monitoringTerminationAction()

    List(describeAction, monitoringAction, terminationAction)
  }

  private def getTaskInputs(inputs: Map[String, WomValue]): String = {
    inputs.transform((name, v) => Map(
      "name" -> name,
      "type" -> v.womType.stableName,
      "value" -> (v match {
        case v: WomFile => v.sizeOption map { _.toString } getOrElse "" // array ? directory ?
        case v: WomValue => v.valueString
      }),
    )).values.toJson.toString
  }

  private def getTaskDisks(disks: Seq[PipelinesApiAttachedDisk], mounts: List[Mount]): String = {
    mounts.flatMap {
      mount => disks
        .filter { disk => disk.name == mount.getDisk }
        .map { disk => Map(
          "type" -> disk.diskType.googleTypeName,
          "path" -> mount.getPath,
        )}
    }.toJson.toString
  }

  def monitoringActions(createPipelineParameters: CreatePipelineParameters, mounts: List[Mount]): List[Action] = {
    val workflowOptions = createPipelineParameters.jobDescriptor.workflowDescriptor.workflowOptions

    lazy val config = workflowOptions.get(WorkflowOptionKeys.MonitoringConfig) getOrElse ""

    workflowOptions.get(WorkflowOptionKeys.MonitoringImage).toOption match {
      case Some(image) => monitoringAction(createPipelineParameters, image, config, mounts)
      case None => List.empty
    }
  }
}
