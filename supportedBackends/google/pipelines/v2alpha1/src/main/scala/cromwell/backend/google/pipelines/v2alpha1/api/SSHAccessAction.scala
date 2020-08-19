package cromwell.backend.google.pipelines.v2alpha1.api

import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import cromwell.backend.google.pipelines.common.action.ActionUtils
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters

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
      .setEntrypoint(ActionUtils.sshEntryPoint)
      .setPortMappings(ActionUtils.sshPortMappings.asJava)
      .setFlags(List(ActionFlag.RunInBackground.toString).asJava)
      .setMounts(mounts.asJava)

    val describeAction = ActionBuilder.describeDocker("ssh access action", sshAction)

    List(describeAction, sshAction)
  }

}
