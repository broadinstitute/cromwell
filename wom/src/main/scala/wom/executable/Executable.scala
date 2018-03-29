package wom.executable

import cats.data.NonEmptyList
import cats.data.Validated.Invalid
import cats.syntax.apply._
import cats.syntax.validated._
import common.Checked
import common.validation.ErrorOr._
import shapeless.Coproduct
import wom.callable.ExecutableCallable
import wom.executable.Executable.ResolvedExecutableInputs
import wom.executable.ExecutableValidation._
import wom.expression.IoFunctionSet
import wom.graph.Graph.ResolvedExecutableInput
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.types.WomType
import wom.values.WomValue

object Executable {

  /*
    * Function provided by each language, that takes the raw input file as a String and returns a Checked[ParsedInputMap]
    * Each entry of the map is an input found in the file.
    * The key is a string representation of the input. It must be be equal to the name of the matching GraphInputNode.
    * The value is a function which given a WomType, attempts to coerce the input value to that type.
    * Thanks to this level of indirection, the logic that links graph input nodes to input values still resides in WOM,
    * which 1) abstracts it away and 2) guarantees that the linking mechanism is the same regardless of the language.
    * At the same time each language can parse the input file however it wants.
   */
  type InputParsingFunction = String => Checked[ParsedInputMap]
  type ParsedInputMap = Map[String, DelayedCoercionFunction]
  type DelayedCoercionFunction = WomType => ErrorOr[WomValue]
  
  /*
    * Maps output ports from graph input nodes to ResolvedExecutableInput
   */
  type ResolvedExecutableInputs = Map[OutputPort, ResolvedExecutableInput]

  def withInputs(entryPoint: ExecutableCallable, inputParsingFunction: InputParsingFunction, inputFile: Option[String], ioFunctions: IoFunctionSet, strictValidation: Boolean): Checked[Executable] = {
    validateExecutable(entryPoint, inputParsingFunction, parseGraphInputs(strictValidation), inputFile, ioFunctions)
  }

  /**
    * Given the graph and the Map[String, DelayedCoercionFunction], attempts to find a value in the map for each ExternalGraphInputNode of the graph
    */
  private def parseGraphInputs(strictValidation: Boolean)(graph: Graph, inputCoercionMap: Map[String, DelayedCoercionFunction], ioFunctions: IoFunctionSet): ErrorOr[ResolvedExecutableInputs] = {
    def fromInputMapping(gin: ExternalGraphInputNode): Option[ErrorOr[ResolvedExecutableInput]] = {
      inputCoercionMap
        .get(gin.nameInInputSet)
        .map { _.apply(gin.womType)
          .map(gin.valueMapper(ioFunctions)(_))
          .map(Coproduct[ResolvedExecutableInput](_))
        }
    }

    def fallBack(gin: ExternalGraphInputNode): ErrorOr[ResolvedExecutableInput] = gin match {
      case required: RequiredGraphInputNode => s"Required workflow input '${required.nameInInputSet}' not specified".invalidNel
      case optionalWithDefault: OptionalGraphInputNodeWithDefault => Coproduct[ResolvedExecutableInput](optionalWithDefault.default).validNel
      case optional: OptionalGraphInputNode => Coproduct[ResolvedExecutableInput](optional.womType.none: WomValue).validNel
    }

    val providedInputsValidation = graph.inputNodes.collect({
      case gin: ExternalGraphInputNode =>
        // The compiler needs the type ascription for some reason
        (gin.singleOutputPort: OutputPort) -> fromInputMapping(gin).getOrElse(fallBack(gin))
    }).toMap.sequence

    val unwantedInputs = if (strictValidation) inputCoercionMap.keySet.diff(graph.externalInputNodes.map(_.nameInInputSet)) else Set.empty

    val wantedInputsValidation: ErrorOr[Unit] = NonEmptyList.fromList(unwantedInputs.toList) match {
      case None => ().validNel
      case Some(unwanteds) => Invalid(unwanteds.map(unwanted => s"WARNING: Unexpected input provided: $unwanted"))
    }

    (providedInputsValidation, wantedInputsValidation) mapN { (providedInputs, _) => providedInputs }
  }
}

/**
  * Wraps a callable in an executable object.
  * @param entryPoint callable that this executable wraps
  * @param resolvedExecutableInputs Resolved values for the ExternalGraphInputNodes of the entryPoint's graph
  */
final case class Executable(entryPoint: ExecutableCallable, resolvedExecutableInputs: ResolvedExecutableInputs) {
  val graph: Graph = entryPoint.graph
}
