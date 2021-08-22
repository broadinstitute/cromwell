package cromwell.backend.google.pipelines.v2beta.api

import com.google.api.services.lifesciences.v2beta.model.{Action, Mount}
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.v2beta.LifeSciencesFactory

trait CheckpointingAction {
  def checkpointingSetupActions(createPipelineParameters: CreatePipelineParameters,
                                mounts: List[Mount]
                            ): List[Action] =
    createPipelineParameters.runtimeAttributes.checkpointFilename map { checkpointFilename =>
      val checkpointingImage = LifeSciencesFactory.CloudSdkImage
      val checkpointingCommand = createPipelineParameters.checkpointingConfiguration.checkpointingCommand(checkpointFilename, ActionCommands.multiLineBinBashCommand)
      val checkpointingEnvironment = Map.empty[String, String]

      // Initial sync from cloud:
      val initialCheckpointSyncAction = ActionBuilder.cloudSdkShellAction(
        createPipelineParameters.checkpointingConfiguration.localizePreviousCheckpointCommand(checkpointFilename)
      )(mounts = mounts)
      val describeInitialCheckpointingSyncAction = ActionBuilder.describeDocker("initial checkpointing sync", initialCheckpointSyncAction)

      // Background upload action:
      val backgroundCheckpointingAction = ActionBuilder.backgroundAction(
        image = checkpointingImage,
        command = checkpointingCommand,
        environment = checkpointingEnvironment,
        mounts = mounts
      )
      val describeBackgroundCheckpointingAction = ActionBuilder.describeDocker("begin checkpointing background action", backgroundCheckpointingAction)

      List(describeInitialCheckpointingSyncAction, initialCheckpointSyncAction, describeBackgroundCheckpointingAction, backgroundCheckpointingAction)
    } getOrElse(Nil)

  def checkpointingShutdownActions(createPipelineParameters: CreatePipelineParameters): List[Action] = {
    createPipelineParameters.runtimeAttributes.checkpointFilename map { checkpointFilename =>
      val terminationAction = ActionBuilder.terminateBackgroundActionsAction()
      val describeTerminationAction = ActionBuilder.describeDocker("terminate checkpointing action", terminationAction)

      val deleteCheckpointAction = ActionBuilder.gcsFileDeletionAction(createPipelineParameters.checkpointingConfiguration.checkpointFileCloud(checkpointFilename))
      val deleteTmpCheckpointAction = ActionBuilder.gcsFileDeletionAction(createPipelineParameters.checkpointingConfiguration.tmpCheckpointFileCloud(checkpointFilename))

      List(describeTerminationAction, terminationAction, deleteCheckpointAction, deleteTmpCheckpointAction)
    } getOrElse(Nil)
  }
}
