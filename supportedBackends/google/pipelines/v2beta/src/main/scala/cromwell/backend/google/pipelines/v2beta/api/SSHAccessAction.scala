package cromwell.backend.google.pipelines.v2beta.api

import com.google.api.services.lifesciences.v2beta.model.{Action, Mount}

import cromwell.backend.google.pipelines.common.WorkflowOptionKeys
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters

import scala.collection.JavaConverters._

trait SSHAccessAction {

  def sshAccessActions(createPipelineParameters: CreatePipelineParameters, mounts: List[Mount]) : List[Action] = {
    val workflowOptions = createPipelineParameters.jobDescriptor.workflowDescriptor.workflowOptions

    workflowOptions.getBoolean(WorkflowOptionKeys.EnableSSHAccess).toOption match {
      case Some(true) => sshAccessAction(mounts)
      case _ => List.empty
    }
  }

  private def sshAccessAction(mounts: List[Mount]): List[Action] = {
    val ports = Map("22" -> new Integer(22))

    val sshAction = ActionBuilder.withImage("gcr.io/cloud-genomics-pipelines/tools")
      .setEntrypoint("ssh-server")
      .setRunInBackground(true)
      .setMounts(mounts.asJava)
      .setPortMappings(mapAsJavaMap(ports))

    List(sshAction)

  }

}
