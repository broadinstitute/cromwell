package cromwell

import cromwell.binding.values.WdlValue
import cromwell.engine.WorkflowId

/**
 * ==WDL Bindings for Scala==
 *
 * This package contains an API to convert WDL ASTs into native Scala objects.
 *
 * Clients of this package should be able to interface with a WDL file without ever coming in direct contact with an Abstract Syntax Tree
 *
 * Internally, this package is built on top of [[cromwell.parser]].
 */

package object binding {
  type WdlSource = String
  type WdlJson = String
  type WorkflowRawInputs = Map[FullyQualifiedName, Any]
  type WorkflowCoercedInputs = Map[FullyQualifiedName, WdlValue]
  type WorkflowOutputs = Map[FullyQualifiedName, WdlValue]
  type FullyQualifiedName = String
  type LocallyQualifiedName = String
  type CallInputs = Map[String, WdlValue]
  type CallOutputs = Map[LocallyQualifiedName, WdlValue]
  type HostInputs = Map[String, WdlValue]

  type ImportResolver = String => WdlSource

  /**
   * Core data identifying a workflow including its unique ID, its namespace, and strongly typed inputs.
   */
  case class WorkflowDescriptor(id: WorkflowId, namespace: NamespaceWithWorkflow, wdlSource: WdlSource, wdlJson: WdlJson, actualInputs: WorkflowCoercedInputs) {
    val name = namespace.workflow.name
    val shortId = id.toString.split("-")(0)
  }
}
