package cromwell

import com.typesafe.config.Config
import cromwell.core.{WorkflowId, WorkflowOptions}
import wdl4s._
import wdl4s.values.WdlValue

package object backend {

  trait JobKey {
    def scope: Scope
    def index: Option[Int]
    def attempt: Int
    def tag: String
  }

  /**
    * For uniquely identifying a job which has been or will be sent to the backend.
    */
  case class BackendJobDescriptorKey(call: Call, index: Option[Int], attempt: Int) extends JobKey {
    val scope = call
    val tag = s"${call.fullyQualifiedName}:$index:$attempt"
  }

  /**
    * For passing to a BackendWorkflowActor for job execution or recovery
    */
  case class BackendJobDescriptor(descriptor: BackendWorkflowDescriptor,
                                  key: BackendJobDescriptorKey,
                                  symbolMap: Map[FullyQualifiedName, WdlValue]) {
    val call = key.call
  }

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
