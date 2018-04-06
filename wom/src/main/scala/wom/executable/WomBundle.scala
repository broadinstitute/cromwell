package wom.executable

import common.Checked
import common.validation.Checked._
import wom.callable.{Callable, CallableTaskDefinition, ExecutableCallable, WorkflowDefinition}
import wom.types.WomType

/**
  * Represents a set of static WOM items that might be imported from a file.
  */
final case class WomBundle(primaryCallable: Option[Callable],
                           allCallables: Set[Callable],
                           typeAliases: Map[String, WomType]) {
  def toExecutableCallable: Checked[ExecutableCallable] = primaryCallable match {
    case Some(w: WorkflowDefinition) => w.validNelCheck
    case Some(c: CallableTaskDefinition) => c.toExecutable.toEither
    case Some(other) => s"Cannot convert WOM bundle to executable. Primary callable was an unknown type ${other.getClass.getSimpleName}.".invalidNelCheck
    case None => s"Cannot convert WOM bundle to executable. No primary callable was available.".invalidNelCheck
  }
}
