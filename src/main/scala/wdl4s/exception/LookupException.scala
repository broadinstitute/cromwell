package wdl4s.exception

import lenthall.exception.ThrowableAggregation
import wdl4s.{Declaration, GraphNode}
import wdl4s.types.WdlType

sealed trait LookupException { this: Exception => }

final case class OptionalNotSuppliedException(operationName: String) extends Exception(s"Sorry! Operation $operationName is not supported on empty optional values. You might resolve this using select_first([optional, default]) to guarantee that you have a filled value.")

/**
  * When a variable does not reference any known scope in the namespace.
  */
sealed abstract case class VariableNotFoundException(variable: String, quoteName: Boolean = true) extends Exception(s"Variable $variable not found") with LookupException

object VariableNotFoundException {
  def apply(variable: String): VariableNotFoundException = new VariableNotFoundException(s"'$variable'") {}
  def apply(variable: String, variableType: WdlType): VariableNotFoundException= new VariableNotFoundException(s"'$variable': ${variableType.toWdlString}") {}
  def apply(declaration: Declaration): VariableNotFoundException = VariableNotFoundException.apply(declaration.fullyQualifiedName, declaration.wdlType)
}

/**
  * When an unexpected exception occurred while attempting to resolve a variable.
  * Might act a single exception or aggregate multiple exceptions.
  */
class VariableLookupException(override val exceptionContext: String, override val throwables: List[Throwable] = List.empty) extends RuntimeException with ThrowableAggregation with LookupException

/**
  * Raised when attempting to resolve a variable in a scatter but the index could not be found.
  */
final case class ScatterIndexNotFound(message: String) extends VariableLookupException(message)

/**
  * Raised when attempting to resolve an output variable but the output resolver failed to return a value.
  */
final case class OutputVariableLookupException(node: GraphNode, index: Option[Int]) extends VariableLookupException(s"Could not find outputs for call ${node.fullyQualifiedName} at index $index")