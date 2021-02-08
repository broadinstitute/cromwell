package cromwell.backend.google.pipelines.v2beta.api

import com.google.api.services.lifesciences.v2beta.model.{Action, Mount}
import cromwell.backend.google.pipelines.common.action.ActionUtils
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.v2beta.api.ActionBuilder._

import scala.collection.JavaConverters._

trait SSHAccessAction {

  def sshAccessActions(createPipelineParameters: CreatePipelineParameters, mounts: List[Mount]) : List[Action] = {
    if (createPipelineParameters.enableSshAccess) {
      sshAccessAction(mounts)
    } else {
      Nil
    }
  }

  private def sshAccessAction(mounts: List[Mount]): List[Action] = {
    val sshAction = ActionBuilder.withImage(ActionUtils.sshImage)
      .withEntrypointCommand(ActionUtils.sshEntryPoint)
      .setPortMappings(ActionUtils.sshPortMappings.asJava)
      .setRunInBackground(true)
      .setMounts(mounts.asJava)

    val describeAction = ActionBuilder.describeDocker("ssh access action", sshAction)

    List(describeAction, sshAction)
  }

}
