package cromwell.binding

import cromwell.binding.Workflow.FQName
import cromwell.parser.WdlParser

object Workflow {
  type FQName = String
}

// Type hierarchy for WDLType with overrides of isCompatible?
trait WDLType {
  def isCompatible(value: Any): Boolean
  def checkCompatible(value: Any) = {
    if (!isCompatible(value)) throw new UnsupportedOperationException
  }
}

class WDLExpression(val ast: WdlParser.Ast) {
  def evaluate: WDLValue = ???
}

class WDLValue(val value: Any, val wdlType: WDLType) {
  wdlType.checkCompatible(value)
}

class CommandPart
class StringCommandPart(literal: String) extends CommandPart
class ParameterCommandPart(wdlType: WDLType, name: String) extends CommandPart

class Command {
  def parts: Seq[CommandPart] = ???
  def inputs: Map[String, WDLType] = ???
  def instantiate(parameters: Map[String, WDLValue]): String = ???
}

class TaskOutput(val name: String, wdlType: WDLType, expression: WDLExpression) {
  def evaluate: WDLValue = {
    val evaluation = expression.evaluate
    wdlType.checkCompatible(evaluation)
    evaluation
  }
}

class Task {
  def name: String = ???
  def command: Command = ???
  def inputs: Map[String, WDLType] = command.inputs
  def outputs: Set[TaskOutput] = ???
}

class Call {
  def name: String = ???
  def task: Task = ???

  /** Mapping from the input section of the call to FQNs of outputs of
    * previous calls. */
  def inputMappings: Map[String, FQName] = ???
  /**
   * Mapping from fully-qualified names to the WDLTypes which a caller must specify.
   */
  def unsatisfiedInputs: Map[String, WDLType] = ???
}


class Workflow {
  def name: String = ???
  
  def calls: Seq[Call] = ???

  def inputs: Map[FQName, WDLType] = ???
  
  def outputs: Map[FQName, WDLType] = ???

}
