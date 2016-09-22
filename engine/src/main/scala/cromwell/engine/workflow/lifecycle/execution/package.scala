package cromwell.engine.workflow.lifecycle

import wdl4s._

package object execution {
  def splitFqn(fullyQualifiedName: FullyQualifiedName): (String, String) = {
    val lastIndex = fullyQualifiedName.lastIndexOf(".")
    (fullyQualifiedName.substring(0, lastIndex), fullyQualifiedName.substring(lastIndex + 1))
  }
}
