package cromwell.backend.google.pipelines.batch.runnable

import cromwell.backend.google.pipelines.batch.GcpBatchRequestFactory.CreatePipelineParameters

import com.google.cloud.batch.v1.{Runnable, Volume}

trait UserRunnable {

  // add in mounts?
  def userRunnables(createPipelineParameters: CreatePipelineParameters, volumes: List[Volume]): List[Runnable] = {
    val userRunnable = RunnableBuilder.userRunnable(
      createPipelineParameters.dockerImage,
      createPipelineParameters.commandScriptContainerPath.pathAsString,
      volumes,
      createPipelineParameters.jobShell,
      //createPipelineParameters.privateDockerKeyAndEncryptedToken,
      //createPipelineParameters.fuseEnabled
    )

    val describeAction = RunnableBuilder.describeDocker("user action", userRunnable.build)

    List(describeAction, userRunnable.build)
  }

}
