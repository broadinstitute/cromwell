package cromwell

import java.util.UUID

import cromwell.binding.values.WdlValue

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
  type CallInputs = Map[String, WdlValue]
  type CallOutputs = Map[FullyQualifiedName, WdlValue]
  type HostInputs = Map[String, WdlValue]

  type ImportResolver = String => WdlSource

  /**
   * Core data identifying a workflow including its unique ID, its namespace, and strongly typed inputs.
   */
  case class WorkflowDescriptor(id: UUID, namespace: NamespaceWithWorkflow, wdlSource: WdlSource, wdlJson: WdlJson, actualInputs: WorkflowCoercedInputs) {
    val name = namespace.workflow.name
    val shortId = id.toString.split("-")(0)
  }
}
