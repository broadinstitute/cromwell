package cromwell.backend.google.pipelines.v2alpha1.api

import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}

trait MonitoringAction {
  def monitoringActions(mounts: List[Mount]): List[Action] = {
    val monitoringAction = ActionBuilder.monitoringAction(mounts)

    val describeAction = ActionBuilder.describeDocker("monitoring action", monitoringAction)

    List(describeAction, monitoringAction)
  }
}
