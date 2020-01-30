package cromwell.backend.google.pipelines.v2alpha1.api

import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import cromwell.backend.google.pipelines.common.WorkflowOptionKeys
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters

import scala.collection.JavaConverters._

trait RemoteAccessAction {

  def remoteAccessActions(createPipelineParameters: CreatePipelineParameters, mounts: List[Mount]) : List[Action] = {
    val workflowOptions = createPipelineParameters.jobDescriptor.workflowDescriptor.workflowOptions

    workflowOptions.getBoolean(WorkflowOptionKeys.EnableRemoteAccess).toOption match {
      case Some(true) => remoteAccessAction(mounts)
      case _ => List.empty
    }
  }

  private def remoteAccessAction(mounts: List[Mount]): List[Action] = {
    val ports = Map("22" -> new Integer(22))

    val sshAction = ActionBuilder.withImage("gcr.io/cloud-genomics-pipelines/tools")
      .setEntrypoint("ssh-server")
      .setFlags(List(ActionFlag.RunInBackground.toString).asJava)
      .setMounts(mounts.asJava)
      .setPortMappings(mapAsJavaMap(ports))

    List(sshAction)

  }

}
