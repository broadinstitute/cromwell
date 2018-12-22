package cromwell.backend.google.pipelines.v2alpha1.api

import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters

trait MonitoringAction {
  object Env {
    /**
      * Name of an environmental variable
      */
    val WorkflowId = "WORKFLOW_ID"
    val TaskCallName = "TASK_CALL_NAME"
    val TaskCallIndex = "TASK_CALL_INDEX"
    val TaskCallAttempt = "TASK_CALL_ATTEMPT"
    val DiskMounts = "DISK_MOUNTS"
  }

  def monitoringActions(createPipelineParameters: CreatePipelineParameters, mounts: List[Mount]): List[Action] = {
    val job = createPipelineParameters.jobDescriptor

    val environment = Map(
      Env.WorkflowId -> job.workflowDescriptor.id.toString,
      Env.TaskCallName -> job.taskCall.localName,
      Env.TaskCallIndex -> (job.key.index map { _.toString } getOrElse "NA"),
      Env.TaskCallAttempt -> job.key.attempt.toString,
      Env.DiskMounts -> mounts.map(_.getPath).mkString(" "),
    )

    val monitoringAction = ActionBuilder.monitoringAction(environment, mounts)

    val describeAction = ActionBuilder.describeDocker("monitoring action", monitoringAction)

    List(describeAction, monitoringAction)
  }
}
