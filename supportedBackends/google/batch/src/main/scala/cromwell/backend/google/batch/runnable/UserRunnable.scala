package cromwell.backend.google.batch.runnable

import com.google.cloud.batch.v1.{Runnable, Volume}
import cromwell.backend.google.batch.api.GcpBatchRequestFactory.CreatePipelineParameters


trait UserRunnable {

  def userRunnables(createPipelineParameters: CreatePipelineParameters, volumes: List[Volume]): List[Runnable] = {

    println(f"job shell ${createPipelineParameters.jobShell}")
    println(f"script container path ${createPipelineParameters.commandScriptContainerPath}")

    val userRunnable = RunnableBuilder.userRunnable(
      docker = createPipelineParameters.dockerImage,
      scriptContainerPath = createPipelineParameters.commandScriptContainerPath.pathAsString,
      // TODO: determine path for shell and scatter
      //jobShell = createPipelineParameters.jobShell,
      jobShell = "/bin/bash",
      volumes = volumes,
      dockerhubCredentials = createPipelineParameters.dockerhubCredentials
      // not necessary for now
      //createPipelineParameters.privateDockerKeyAndEncryptedToken,
      //createPipelineParameters.fuseEnabled
    )

    val describeRunnable = RunnableBuilder.describeDocker("user runnable", userRunnable)
    List(describeRunnable, userRunnable).map(_.build)
  }
}
