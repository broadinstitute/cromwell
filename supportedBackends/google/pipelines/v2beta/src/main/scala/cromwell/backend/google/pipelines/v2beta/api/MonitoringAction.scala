package cromwell.backend.google.pipelines.v2beta.api

import com.google.api.services.lifesciences.v2beta.model.{Action, Mount}
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.{CreatePipelineParameters, MountsToEnv}
import cromwell.backend.google.pipelines.v2beta.api.ActionBuilder.Labels.{Key, Value}

trait MonitoringAction {
  private def monitoringAction(job: BackendJobDescriptor,
                               monitoringImage: String,
                               monitoringImageCommand: List[String],
                               monitoringImageEnvironment: MountsToEnv,
                               mounts: List[Mount],
                              ): List[Action] = {
    val monitoringAction = ActionBuilder.monitoringAction(
      monitoringImage,
      monitoringImageCommand,
      monitoringImageEnvironment(mounts.map(_.getPath)),
      mounts
    )

    val describeAction = ActionBuilder.describeDocker("monitoring action", monitoringAction).setRunInBackground(true)

    val terminationAction = ActionBuilder.monitoringTerminationAction()

    val describeTerminationAction = ActionBuilder.describeDocker("terminate monitoring action", monitoringAction)

    List(describeAction, monitoringAction, describeTerminationAction, terminationAction)
  }

  def monitoringActions(createPipelineParameters: CreatePipelineParameters, mounts: List[Mount]): List[Action] = {
    createPipelineParameters.monitoringImage match {
      case Some(image) =>
        monitoringAction(
          job = createPipelineParameters.jobDescriptor,
          monitoringImage = image,
          monitoringImageCommand = createPipelineParameters.monitoringImageCommand,
          monitoringImageEnvironment = createPipelineParameters.monitoringImageEnvironment,
          mounts = mounts,
        )
      case None => List.empty
    }
  }

  def monitoringPreambleActions(createPipelineParameters: CreatePipelineParameters,
                                mounts: List[Mount]
                               )(implicit gcsTransferConfiguration: GcsTransferConfiguration): List[Action] = {
    createPipelineParameters.monitoringImageScript match {
      case Some(script) =>
        val command = ActionCommands.localizeFile(script, createPipelineParameters.monitoringImageScriptContainerPath)
        val labels = Map(Key.Tag -> Value.Localization)
        val localizeScriptAction = ActionBuilder.cloudSdkShellAction(command)(mounts = mounts, labels = labels)
        val describeLocalizeScriptAction = ActionBuilder.describeDocker(
          "localizing monitoring image script action",
          localizeScriptAction,
        )
        List(describeLocalizeScriptAction, localizeScriptAction)
      case None => List.empty
    }
  }
}
