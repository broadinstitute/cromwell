package cwl

import cats.data.Validated.{Invalid, Valid}
import cats.data.{NonEmptyList, OptionT}
import cats.syntax.traverse._
import cats.syntax.validated._
import common.Checked
import common.validation.ErrorOr._
import cwl.CommandLineTool._
import cwl.CwlVersion._
import cwl.command.ParentName
import io.circe.Json
import shapeless.syntax.singleton._
import shapeless.{:+:, CNil, Coproduct, Poly1, Witness}
import wom.callable.CommandTaskDefinition.{EvaluatedOutputs, OutputFunctionResponse}
import wom.callable.{Callable, CallableTaskDefinition}
import wom.expression.{IoFunctionSet, ValueAsAnExpression, WomExpression}
import wom.graph.GraphNodePort.OutputPort
import wom.types.{WomOptionalType, WomStringType}
import wom.values.{WomEvaluatedCallInputs, WomGlobFile, WomString, WomValue}
import wom.{CommandPart, RuntimeAttributes}

import scala.concurrent.ExecutionContext
import scala.math.Ordering

/**
  * @param `class` This _should_ always be "CommandLineTool," however the spec does not -er- specify this.
  */
case class CommandLineTool private(
                                    inputs: Array[CommandInputParameter],
                                    outputs: Array[CommandOutputParameter],
                                    `class`: Witness.`"CommandLineTool"`.T,
                                    id: String,
                                    requirements: Option[Array[Requirement]],
                                    hints: Option[Array[Hint]],
                                    label: Option[String],
                                    doc: Option[String],
                                    cwlVersion: Option[CwlVersion],
                                    baseCommand: Option[BaseCommand],
                                    arguments: Option[Array[CommandLineTool.Argument]],
                                    stdin: Option[StringOrExpression],
                                    stderr: Option[StringOrExpression],
                                    stdout: Option[StringOrExpression],
                                    successCodes: Option[Array[Int]],
                                    temporaryFailCodes: Option[Array[Int]],
                                    permanentFailCodes: Option[Array[Int]]) extends Tool {

  def asCwl: Cwl = Coproduct[Cwl](this)

  /*
   * The command template is built following the rules described here: http://www.commonwl.org/v1.0/CommandLineTool.html#Input_binding
   * - The baseCommand goes first
   * - Then the arguments are assigned a sorting key and transformed into a CommandPart
   * - Finally the inputs are folded one by one into a CommandPartsList
   * - arguments and inputs CommandParts are sorted according to their sort key
   */
  private [cwl] def buildCommandTemplate(expressionLib: ExpressionLib)(inputValues: WomEvaluatedCallInputs): ErrorOr[List[CommandPart]] = {
    import cats.instances.list._
    import cats.syntax.traverse._

    val baseCommandPart = baseCommand.toList.flatMap(_.fold(BaseCommandToCommandParts))

    val argumentsParts: CommandPartsList =
    // arguments is an Option[Array[Argument]], the toList.flatten gives a List[Argument]
      arguments.toList.flatten
        // zip the index because we need it in the sorting key
        .zipWithIndex.foldLeft(CommandPartsList.empty)({
        case (commandPartsList, (argument, index)) =>
          val part = argument.fold(ArgumentToCommandPart).apply(expressionLib)
          // Get the position from the binding if there is one
          val position = argument.select[ArgumentCommandLineBinding].flatMap(_.position)
            .map(Coproduct[StringOrInt](_)).getOrElse(DefaultPosition)

          // The key consists of the position followed by the index
          val sortingKey = CommandBindingSortingKey(List(position, Coproduct[StringOrInt](index)))

          commandPartsList :+ SortKeyAndCommandPart(sortingKey, part)
      })

    val inputBindingsCommandParts: ErrorOr[List[SortKeyAndCommandPart]] = inputs.toList.flatTraverse[ErrorOr, SortKeyAndCommandPart]({
      import cats.syntax.validated._

      inputParameter =>
        val parsedName = FullyQualifiedName(inputParameter.id)(ParentName.empty).id
        val womType = inputParameter.`type`.map(_.fold(MyriadInputTypeToWomType)).getOrElse(WomStringType)

        val defaultValue = inputParameter.default.map(_.fold(InputParameter.DefaultToWomValuePoly).apply(womType))

        inputValues
          .collectFirst({ case (inputDefinition, womValue) if inputDefinition.name == parsedName => womValue.validNel })
          .orElse(defaultValue) match {
          case Some(Valid(value)) =>
            // See http://www.commonwl.org/v1.0/CommandLineTool.html#Input_binding
            lazy val initialKey = CommandBindingSortingKey.empty
              .append(inputParameter.inputBinding, Coproduct[StringOrInt](parsedName))

            inputParameter.`type`.toList.flatMap(_.fold(MyriadInputTypeToSortedCommandParts).apply(inputParameter.inputBinding, value, initialKey.asNewKey, expressionLib)).validNel
          case Some(Invalid(errors)) => Invalid(errors)
          case None => s"Could not find an input value for input $parsedName in ${inputValues.prettyString}".invalidNel
        }
    })

    inputBindingsCommandParts map { parts =>
      baseCommandPart ++ (argumentsParts ++ parts).sorted.map(_.commandPart)
    }
  }

  // This seems like it makes sense only for CommandLineTool and hence is not abstracted in Tool. If this assumption is wrong it could be moved up.
  private def environmentDefs(requirementsAndHints: List[Requirement], expressionLib: ExpressionLib): ErrorOr[Map[String, WomExpression]] = {
    // For environment variables we need to make sure that we aren't being asked to evaluate expressions from a containing
    // workflow step or its containing workflow or anything containing the workflow. The current structure of this code
    // is not prepared to evaluate those expressions. Actually this is true for attributes too and we're totally not
    // checking for this condition there. Blurgh.
    // TODO CWL: for runtime attributes, detect unevaluatable expressions in the containment hierarchy.

    // This traverses all `EnvironmentDef`s within all `EnvVarRequirement`s. The spec doesn't appear to say how to handle
    // duplicate `envName` keys in a single array of `EnvironmentDef`s; this code gives precedence to the last occurrence.
    val allEnvVarDefs = for {
      req <- requirementsAndHints
      envVarReq <- req.select[EnvVarRequirement].toList
      // Reverse the defs within an env var requirement so that when we fold from the right below the later defs
      // will take precedence over the earlier defs.
      envDef <- envVarReq.envDef.toList.reverse
    } yield envDef

    // Compact the `EnvironmentDef`s. Don't convert to `WomExpression`s yet, the `StringOrExpression`s need to be
    // compared to the `EnvVarRequirement`s that were defined on this tool.
    val effectiveEnvironmentDefs = allEnvVarDefs.foldRight(Map.empty[String, StringOrExpression]) {
      case (envVarReq, envVarMap) => envVarMap + (envVarReq.envName -> envVarReq.envValue)
    }

    // These are the effective environment defs irrespective of where they were found in the
    // Run / WorkflowStep / Workflow containment hierarchy.
    val effectiveExpressionEnvironmentDefs = effectiveEnvironmentDefs filter { case (_, expr) => expr.select[Expression].isDefined }

    // These are only the environment defs defined on this tool.
    val cltRequirements = requirements.toList.flatten ++ hints.toList.flatten.flatMap(_.select[Requirement])
    val cltEnvironmentDefExpressions = (for {
      cltEnvVarRequirement <- cltRequirements flatMap {
        _.select[EnvVarRequirement]
      }
      cltEnvironmentDef <- cltEnvVarRequirement.envDef.toList
      expr <- cltEnvironmentDef.envValue.select[Expression].toList
    } yield expr).toSet

    // If there is an expression in an effective environment def that wasn't defined on this tool then error out since
    // there isn't currently a way of evaluating it.
    val unevaluatableEnvironmentDefs = for {
      (name, stringOrExpression) <- effectiveExpressionEnvironmentDefs.toList
      expression <- stringOrExpression.select[Expression].toList
      if !cltEnvironmentDefExpressions.contains(expression)
    } yield name

    unevaluatableEnvironmentDefs match {
      case Nil =>
        // No unevaluatable environment defs => keep on truckin'
        effectiveEnvironmentDefs.foldRight(Map.empty[String, WomExpression]) { case ((envName, envValue), acc) =>
          acc + (envName -> envValue.fold(StringOrExpressionToWomExpression).apply(inputNames, expressionLib))
        }.validNel
      case xs =>
        s"Could not evaluate environment variable expressions defined in the call hierarchy of tool $id: ${xs.mkString(", ")}.".invalidNel
    }
  }

  /*
    * Custom evaluation of the outputs of the command line tool (custom as in bypasses the engine provided default implementation).
    * This is needed because the output of a CWL tool might be defined by the presence of a cwl.output.json file in the output directory.
    * In that case, the content of the file should be used as output. This means that otherwise provided output expressions should not be evaluated.
    * This method checks for the existence of this file, and optionally return a Map[OutputPort, WomValue] if the file was found.
    * It ensures all the output ports of the corresponding WomCallNode get a WomValue which is needed for the engine to keep running the workflow properly.
    * If the json file happens to be missing values for one or more output ports, it's a failure.
    * If the file is not present, an empty value is returned and the engine will execute its own output evaluation logic.
    *
    * TODO: use IO instead of Future ?
   */
  final private def outputEvaluationJsonFunction(outputPorts: Set[OutputPort],
                                                 inputs: Map[String, WomValue],
                                                 ioFunctionSet: IoFunctionSet,
                                                 executionContext: ExecutionContext): OutputFunctionResponse = {
    implicit val ec = executionContext
    import cats.instances.future._
    import cats.syntax.either._

    // Convert the parsed json to Wom-usable output values
    def jsonToOutputs(json: Map[String, Json]): Checked[List[(OutputPort, WomValue)]] = {
      import cats.instances.list._

      outputPorts.toList.traverse[ErrorOr, (OutputPort, WomValue)]({ outputPort =>
        // If the type is optional, then we can set the value to none if there's nothing in the json
        def emptyValue = outputPort.womType match {
          case optional: WomOptionalType => Option((outputPort -> optional.none).validNel)
          case _ => None
        }

        json.get(outputPort.name)
          .map(_.foldWith(CwlJsonToDelayedCoercionFunction).apply(outputPort.womType).map(outputPort -> _))
          .orElse(emptyValue)
          .getOrElse(s"Cannot find a value for output ${outputPort.name} in output json $json".invalidNel)
      }).toEither
    }

    // Parse content as json and return output values for each output port
    def parseContent(content: String): EvaluatedOutputs = {
      for {
        parsed <- io.circe.yaml.parser.parse(content).flatMap(_.as[Map[String, Json]]).leftMap(error => NonEmptyList.one(error.getMessage))
        jobOutputsMap <- jsonToOutputs(parsed)
      } yield jobOutputsMap.toMap
    }

    for {
      // Glob for "cwl.output.json"
      outputJsonGlobs <- OptionT.liftF { ioFunctionSet.glob(CwlOutputJson) }
      // There can only be 0 or 1, so try to take the head of the list
      outputJsonFile <- OptionT.fromOption { outputJsonGlobs.headOption }
      // Read the content using the ioFunctionSet.readFile function
      content <- OptionT.liftF { ioFunctionSet.readFile(outputJsonFile, None, failOnOverflow = false) }
      // parse the content and validate it
      outputs = parseContent(content)
    } yield outputs
  }

  def buildTaskDefinition(taskName: String,
                          inputDefinitions: List[_ <: Callable.InputDefinition],
                          outputDefinitions: List[Callable.OutputDefinition],
                          runtimeAttributes: RuntimeAttributes,
                          requirementsAndHints: List[cwl.Requirement],
                          expressionLib: ExpressionLib): ErrorOr[CallableTaskDefinition] = {

    def redirect(soe: Option[StringOrExpression]): Option[WomExpression] = soe map {
      _.fold(StringOrExpressionToWomExpression).apply(inputNames, expressionLib)
    }

    // This seems like it makes sense only for CommandLineTool and hence is not abstracted in Tool. If this assumption is wrong it could be moved up.
    val adHocFileCreations: Set[WomExpression] = (for {
      requirements <- requirements.getOrElse(Array.empty[Requirement])
      initialWorkDirRequirement <- requirements.select[InitialWorkDirRequirement].toArray
      listing <- initialWorkDirRequirement.listings
    } yield InitialWorkDirFileGeneratorExpression(listing, expressionLib)).toSet[WomExpression]

    environmentDefs(requirementsAndHints, expressionLib) map { environmentExpressions =>
      CallableTaskDefinition(
        taskName,
        buildCommandTemplate(expressionLib),
        runtimeAttributes,
        Map.empty,
        Map.empty,
        outputDefinitions,
        inputDefinitions,
        adHocFileCreation = adHocFileCreations,
        environmentExpressions = environmentExpressions,
        // TODO: This doesn't work in all cases and it feels clunky anyway - find a way to sort that out
        prefixSeparator = "#",
        commandPartSeparator = " ",
        stdinRedirection = redirect(stdin),
        stdoutOverride = redirect(stdout),
        stderrOverride = redirect(stderr),
        additionalGlob = Option(WomGlobFile(CwlOutputJson)),
        customizedOutputEvaluation = outputEvaluationJsonFunction,
        homeOverride = Option(_.outputPath)
      )
    }
  }
}

object CommandLineTool {
  val CwlOutputJson = "cwl.output.json"

  private val DefaultPosition = Coproduct[StringOrInt](0)
  // Elements of the sorting key can be either Strings or Ints
  type StringOrInt = String :+: Int :+: CNil
  object StringOrInt {
    object String { def unapply(stringOrInt: StringOrInt): Option[String] = stringOrInt.select[String] }
    object Int { def unapply(stringOrInt: StringOrInt): Option[Int] = stringOrInt.select[Int] }
  }

  /*
   * The algorithm described here http://www.commonwl.org/v1.0/CommandLineTool.html#Input_binding to sort the command line parts
   * uses a sorting key assigned to each binding. This class represents such a key
   */
  object CommandBindingSortingKey {
    def empty = CommandBindingSortingKey(List.empty, List.empty)
  }
  case class CommandBindingSortingKey(head: List[StringOrInt],
                                      tail: List[StringOrInt] = List.empty) {
    val value = head ++ tail

    def append(binding: Option[CommandLineBinding], name: StringOrInt): CommandBindingSortingKey = binding match {
      // If there's an input binding, add the position to the key (or 0)
      case Some(b) =>
        // The spec is inconsistent about this as it says "If position is not specified, it is not added to the sorting key"
        // but also that the position defaults to 0. cwltool uses 0 when there's no position so we'll do that too.
        val position = b.position.map(Coproduct[StringOrInt](_)) getOrElse DefaultPosition
        copy(head = head :+ position, tail = tail :+ name)
      // Otherwise do nothing
      case None => this
    }

    /**
      * Creates a new key with head and tail combined into the new head.
      */
    def asNewKey = CommandBindingSortingKey(value)
  }

  // Maps a sorting key to its binding
  case class SortKeyAndCommandPart(sortingKey: CommandBindingSortingKey, commandPart: CommandPart)

  type CommandPartsList = List[SortKeyAndCommandPart]

  object CommandPartsList {
    def empty: CommandPartsList = List.empty[SortKeyAndCommandPart]
  }

  // Ordering for CommandBindingSortingKeyElement
  implicit val SortingKeyTypeOrdering: Ordering[StringOrInt] = Ordering.fromLessThan[StringOrInt]({
    // String comparison
    case (StringOrInt.String(s1), StringOrInt.String(s2)) => s1 < s2
    // Int comparison
    case (StringOrInt.Int(i1), StringOrInt.Int(i2)) => i1 < i2
    // Int < String (from the spec: "Numeric entries sort before strings.")
    case (StringOrInt.Int(_), StringOrInt.String(_)) => true
    // String > Int
    case (StringOrInt.String(_), StringOrInt.Int(_)) => false
  })

  // Ordering for a CommandBindingSortingKey
  implicit val SortingKeyOrdering: Ordering[CommandBindingSortingKey] = Ordering.by(_.value.toIterable)

  // Ordering for a CommandPartSortMapping: order by sorting key
  implicit val SortKeyAndCommandPartOrdering: Ordering[SortKeyAndCommandPart] = Ordering.by(_.sortingKey)

  def apply(inputs: Array[CommandInputParameter] = Array.empty,
            outputs: Array[CommandOutputParameter] = Array.empty,
            id: String,
            requirements: Option[Array[Requirement]] = None,
            hints: Option[Array[Hint]] = None,
            label: Option[String] = None,
            doc: Option[String] = None,
            cwlVersion: Option[CwlVersion] = Option(CwlVersion.Version1),
            baseCommand: Option[BaseCommand] = None,
            arguments: Option[Array[CommandLineTool.Argument]] = None,
            stdin: Option[StringOrExpression] = None,
            stderr: Option[StringOrExpression] = None,
            stdout: Option[StringOrExpression] = None,
            successCodes: Option[Array[Int]] = None,
            temporaryFailCodes: Option[Array[Int]] = None,
            permanentFailCodes: Option[Array[Int]] = None): CommandLineTool =
    CommandLineTool(inputs, outputs, "CommandLineTool".narrow, id, requirements, hints, label, doc, cwlVersion, baseCommand, arguments, stdin, stderr, stdout, successCodes, temporaryFailCodes, permanentFailCodes)

  type BaseCommand = SingleOrArrayOfStrings

  type Argument = ArgumentCommandLineBinding :+: StringOrExpression :+: CNil

  case class CommandInputParameter(
                                    id: String,
                                    label: Option[String] = None,
                                    secondaryFiles: Option[SecondaryFiles] = None,
                                    format: Option[InputParameterFormat] = None, //only valid when type: File
                                    streamable: Option[Boolean] = None, //only valid when type: File
                                    doc: Option[Doc] = None,
                                    inputBinding: Option[InputCommandLineBinding] = None,
                                    default: Option[CwlAny] = None,
                                    `type`: Option[MyriadInputType] = None) extends InputParameter

  case class CommandOutputParameter(
                                     id: String,
                                     label: Option[String] = None,
                                     secondaryFiles: Option[SecondaryFiles] = None,
                                     format: Option[OutputParameterFormat] = None, //only valid when type: File
                                     streamable: Option[Boolean] = None, //only valid when type: File
                                     doc: Option[Doc] = None,
                                     outputBinding: Option[CommandOutputBinding] = None,
                                     `type`: Option[MyriadOutputType] = None) extends OutputParameter
}

object StringOrExpressionToWomExpression extends Poly1 {
  implicit def string: Case.Aux[String, (Set[String], ExpressionLib) => WomExpression] = at[String] { s => (inputNames, expressionLib) =>
    ValueAsAnExpression(WomString(s))
  }

  implicit def expression: Case.Aux[Expression, (Set[String], ExpressionLib) => WomExpression] = at[Expression] { e => (inputNames, expressionLib) =>
    cwl.ECMAScriptWomExpression(e, inputNames, expressionLib)
  }
}
