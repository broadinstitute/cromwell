package cromwell.backend.google.pipelines.v2alpha1.api

import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters

trait MonitoringAction {
  def monitoringSetupActions(createPipelineParameters: CreatePipelineParameters,
                             mounts: List[Mount]
                            )(implicit gcsTransferConfiguration: GcsTransferConfiguration): List[Action] = {

    val monitoringImageScriptActions =
      createPipelineParameters.monitoringImageScript match {
        case Some(script) =>
          val localizeScriptAction =
            ActionBuilder.monitoringImageScriptAction(
              script,
              createPipelineParameters.monitoringImageScriptContainerPath,
              mounts,
            )
          val describeLocalizeScriptAction =
            ActionBuilder.describeDocker(
              "localizing monitoring image script action",
              localizeScriptAction,
            )
          List(describeLocalizeScriptAction, localizeScriptAction)
        case None => Nil
      }

    val monitoringImageActions =
      createPipelineParameters.monitoringImage match {
      case Some(image) =>

        val monitoringImage = image
        val monitoringImageCommand = createPipelineParameters.monitoringImageCommand
        val monitoringImageEnvironment = createPipelineParameters.monitoringImageEnvironment

        val monitoringAction = ActionBuilder.monitoringAction(
          monitoringImage,
          monitoringImageCommand,
          monitoringImageEnvironment(mounts.map(_.getPath)),
          mounts,
        )
        val describeMonitoringAction = ActionBuilder.describeDocker("monitoring action", monitoringAction)
        List(describeMonitoringAction, monitoringAction)

      case None => Nil
    }

    monitoringImageScriptActions ++ monitoringImageActions
  }

  def monitoringShutdownActions(createPipelineParameters: CreatePipelineParameters): List[Action] = {
    createPipelineParameters.monitoringImage match {
      case Some(_) =>
        val terminationAction = ActionBuilder.monitoringTerminationAction()

        val describeTerminationAction = ActionBuilder.describeDocker("terminate monitoring action", terminationAction)

        List(describeTerminationAction, terminationAction)
      case None => Nil
    }
  }
}
