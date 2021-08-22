package cromwell.backend.google.pipelines.v2beta.api

import com.google.api.services.lifesciences.v2beta.model.{Action, Mount}
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters

trait MonitoringAction {
  def monitoringSetupActions(createPipelineParameters: CreatePipelineParameters,
                             mounts: List[Mount]
                            )(implicit gcsTransferConfiguration: GcsTransferConfiguration): List[Action] = {

    val monitoringImageScriptActions =
      createPipelineParameters.monitoringImage.monitoringImageScriptOption match {
        case Some(script) =>
          val localizeScriptAction =
            ActionBuilder.monitoringImageScriptAction(
              script,
              createPipelineParameters.monitoringImage.monitoringImageScriptContainerPath,
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
      createPipelineParameters.monitoringImage.monitoringImageOption match {
        case Some(image) =>

          val monitoringImage = image
          val monitoringImageCommand = createPipelineParameters.monitoringImage.monitoringImageCommand
          val monitoringImageEnvironment = createPipelineParameters.monitoringImage.monitoringImageEnvironment

          val monitoringAction = ActionBuilder.backgroundAction(
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
    createPipelineParameters.monitoringImage.monitoringImageOption match {
      case Some(_) =>
        val terminationAction = ActionBuilder.terminateBackgroundActionsAction()

        val describeTerminationAction = ActionBuilder.describeDocker("terminate monitoring action", terminationAction)

        List(describeTerminationAction, terminationAction)
      case None => Nil
    }
  }
}
