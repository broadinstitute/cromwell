package cromwell

import com.typesafe.config.Config
import cromwell.core.{WorkflowId, WorkflowOptions}
import wdl4s._
import wdl4s.values.WdlValue

package object backend {

  /**
    * For uniquely identifying a job which has been or will be sent to the backend.
    */
  case class BackendJobDescriptorKey(callFqn: String, index: Option[Int], attempt: Int)

  /**
    * For passing to a BackendWorkflowActor for job execution or recovery
    */
  //case class BackendJobDescriptor(key: BackendJobDescriptorKey, call: Call, locallyQualifiedInputs: CallInputs)
  case class BackendJobDescriptor(key: BackendJobDescriptorKey, call: Call, symbolMap: Map[FullyQualifiedName, WdlValue])

  /**
    * For passing to a BackendActor construction time
    */
  case class BackendWorkflowDescriptor(id: WorkflowId,
                                       workflowNamespace: NamespaceWithWorkflow,
                                       inputs: Map[String, WdlValue],
                                       workflowOptions: WorkflowOptions)

  /**
    * For passing to a BackendActor construction time
    */
  case class BackendConfigurationDescriptor(configPath: String, config: Config)
}
