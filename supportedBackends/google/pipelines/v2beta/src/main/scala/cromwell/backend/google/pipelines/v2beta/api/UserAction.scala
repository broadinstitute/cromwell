package cromwell.backend.google.pipelines.v2beta.api

import com.google.api.services.lifesciences.v2beta.model.{Action, Mount}
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters

trait UserAction {
  def userActions(createPipelineParameters: CreatePipelineParameters, mounts: List[Mount]): List[Action] = {
    val userAction = ActionBuilder.userAction(
      createPipelineParameters.dockerImage,
      createPipelineParameters.commandScriptContainerPath.pathAsString,
      mounts,
      createPipelineParameters.jobShell,
      createPipelineParameters.privateDockerKeyAndEncryptedToken,
      createPipelineParameters.fuseEnabled
    )

    val describeAction = ActionBuilder.describeDocker("user action", userAction)

    List(describeAction, userAction)
  }
}
