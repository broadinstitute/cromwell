package cromwell.backend.google.pipelines.batch.runnable

import com.google.cloud.batch.v1.Runnable
import cromwell.backend.google.pipelines.batch.api.GcpBatchRequestFactory.CreatePipelineParameters

trait UserRunnable {

  // add in mounts?
  def userRunnables(createPipelineParameters: CreatePipelineParameters): List[Runnable] = {
    val userRunnable = RunnableBuilder.userRunnable(
      docker = createPipelineParameters.dockerImage,
      // TODO: Alex - This used to be createPipelineParameters.commandScriptContainerPath.pathAsString which is /cromwell_root/script
      // I'm not sure if such a script will include the value from gcpBatchCommand
      command = createPipelineParameters.gcpBatchCommand,
      jobShell = createPipelineParameters.jobShell,
      // not necessary for now
      //createPipelineParameters.privateDockerKeyAndEncryptedToken,
      //createPipelineParameters.fuseEnabled
    )

    val describeAction = RunnableBuilder.describeDocker("user action", userRunnable.build)
    List(describeAction, userRunnable.build)
  }
}
