package wom.executable

import cats.syntax.validated._
import lenthall.Checked
import lenthall.validation.ErrorOr._
import shapeless.Coproduct
import wom.callable.Callable
import wom.executable.Executable.ResolvedExecutableInputs
import wom.executable.ExecutableValidation._
import wom.graph.Graph.ResolvedExecutableInput
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.types.WdlType
import wom.values.WdlValue

object Executable {

  /*
    * Function provided by each language, that takes the raw input file as a String and returns a Checked[ParsedInputMap]
    * Each entry of the map is an input found in the file.
    * The key is a string representation of the input. It must be be equal to the name of the matching GraphInputNode.
    * The value is a function which given a WdlType, attempts to coerce the input value to that type.
    * Thanks to this level of indirection, the logic that links graph input nodes to input values still resides in WOM,
    * which 1) abstracts it away and 2) guarantees that the linking mechanism is the same regardless of the language.
    * At the same time each language can parse the input file however it wants.
   */
  type InputParsingFunction = String => Checked[ParsedInputMap]
  type ParsedInputMap = Map[String, DelayedCoercionFunction]
  type DelayedCoercionFunction = WdlType => ErrorOr[WdlValue]
  
  /*
    * Maps output ports from graph input nodes to ResolvedExecutableInput
   */
  type ResolvedExecutableInputs = Map[OutputPort, ResolvedExecutableInput]

  def withInputs(entryPoint: Callable, inputParsingFunction: InputParsingFunction, inputFile: Option[String]): Checked[Executable] = {
    validateExecutable(entryPoint, inputParsingFunction, parseGraphInputs, inputFile)
  }

  /**
    * Given the graph and the Map[String, DelayedCoercionFunction], attempts to find a value in the map for each ExternalGraphInputNode of the graph
    */
  private def parseGraphInputs(graph: Graph, inputCoercionMap: Map[String, DelayedCoercionFunction]): ErrorOr[ResolvedExecutableInputs] = {
    def fromInputMapping(gin: ExternalGraphInputNode): Option[ErrorOr[ResolvedExecutableInput]] = {
      inputCoercionMap.get(gin.identifier.fullyQualifiedName.value).map(_(gin.womType).map(Coproduct[ResolvedExecutableInput](_)))
    }

    def fallBack(gin: ExternalGraphInputNode): ErrorOr[ResolvedExecutableInput] = gin match {
      case required: RequiredGraphInputNode => s"Required workflow input '${required.identifier.fullyQualifiedName.value}' not specified".invalidNel
      case optionalWithDefault: OptionalGraphInputNodeWithDefault => Coproduct[ResolvedExecutableInput](optionalWithDefault.default).validNel
      case optional: OptionalGraphInputNode => Coproduct[ResolvedExecutableInput](optional.womType.none: WdlValue).validNel
    }

    graph.inputNodes.collect({
      case gin: ExternalGraphInputNode =>
        // The compiler needs the type ascription for some reason
        (gin.singleOutputPort: OutputPort) -> fromInputMapping(gin).getOrElse(fallBack(gin))
    }).toMap.sequence
  }
}

/**
  * Wraps an callable in an executable object. An executable 
  * @param entryPoint callable that this executable wraps
  * @param resolvedExecutableInputs Resolved values for the ExternalGraphInputNodes of the entryPoint's graph
  */
final case class Executable(entryPoint: Callable, resolvedExecutableInputs: ResolvedExecutableInputs) {
  val graph: ErrorOr[Graph] = entryPoint.graph
}
