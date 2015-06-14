package cromwell

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
  type WorkflowRawInputs = Map[FullyQualifiedName, Any]
  type WorkflowCoercedInputs = Map[FullyQualifiedName, WdlValue]
  type WorkflowOutputs = Map[FullyQualifiedName, WdlValue]
  type FullyQualifiedName = String
  type RuntimeAttributes = Map[String, String]
  type RuntimeAttribute = (String, String)

  /**
   * Provides a few convenience methods for specific runtime attribute keys which were defined in the WDL spec
   */
  implicit class EnhancedRuntimeAttributes(val runtimeAttributes: RuntimeAttributes) extends AnyVal {
    def docker: Option[String] = attribute("docker")
    def memory: Option[String] = attribute("memory")
    def serialize: Option[String] = attribute("serialize")
    def attribute(attr: String): Option[String] = runtimeAttributes.get(attr)
  }
}
