package cromwell.backend.google.batch.runnable

import com.google.cloud.batch.v1.{Runnable, Volume}
import cromwell.backend.google.batch.api.GcpBatchRequestFactory.CreateBatchJobParameters

trait CheckpointingRunnable {
  def checkpointingSetupRunnables(createParameters: CreateBatchJobParameters, volumes: List[Volume]): List[Runnable] = {
    val result = createParameters.runtimeAttributes.checkpointFilename
      .map { checkpointFilename =>
        val checkpointingImage = RunnableUtils.CloudSdkImage
        val checkpointingCommand =
          createParameters.checkpointingConfiguration.checkpointingCommand(checkpointFilename,
                                                                           RunnableCommands.multiLineBinBashCommand
          )
        val checkpointingEnvironment = Map.empty[String, String]

        // Initial sync from cloud:
        val initialCheckpointSyncRunnable = RunnableBuilder.cloudSdkShellRunnable(
          createParameters.checkpointingConfiguration.localizePreviousCheckpointCommand(checkpointFilename)
        )(volumes = volumes, flags = List.empty, labels = Map.empty)
        val describeInitialCheckpointingSyncRunnable =
          RunnableBuilder.describeDocker("initial checkpointing sync", initialCheckpointSyncRunnable)

        // Background upload runnable:
        val backgroundCheckpointingRunnable = RunnableBuilder.backgroundRunnable(
          image = checkpointingImage,
          command = checkpointingCommand,
          environment = checkpointingEnvironment,
          volumes = volumes
        )
        val describeBackgroundCheckpointingRunnable =
          RunnableBuilder.describeDocker("begin checkpointing background runnable", backgroundCheckpointingRunnable)

        List(describeInitialCheckpointingSyncRunnable,
             initialCheckpointSyncRunnable,
             describeBackgroundCheckpointingRunnable,
             backgroundCheckpointingRunnable
        )
      }
      .getOrElse(Nil)

    result.map(_.build)
  }

  def checkpointingShutdownRunnables(createParameters: CreateBatchJobParameters,
                                     volumes: List[Volume]
  ): List[Runnable] = {
    val result = createParameters.runtimeAttributes.checkpointFilename
      .map { checkpointFilename =>
        val terminationRunnable = RunnableBuilder.terminateBackgroundRunnablesRunnable()
        val describeTerminationRunnable =
          RunnableBuilder.describeDocker("terminate checkpointing runnable", terminationRunnable)

        val deleteCheckpointRunnable = RunnableBuilder.gcsFileDeletionRunnable(
          createParameters.checkpointingConfiguration.checkpointFileCloud(checkpointFilename),
          volumes
        )
        val deleteTmpCheckpointRunnable = RunnableBuilder.gcsFileDeletionRunnable(
          createParameters.checkpointingConfiguration.tmpCheckpointFileCloud(checkpointFilename),
          volumes
        )

        List(describeTerminationRunnable, terminationRunnable, deleteCheckpointRunnable, deleteTmpCheckpointRunnable)
      }
      .getOrElse(Nil)

    result.map(_.build)
  }
}
