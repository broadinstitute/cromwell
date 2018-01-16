package cwl

import java.nio.file.Paths

import cats.data.Validated.{Invalid, Valid}
import cats.instances.list._
import cats.instances.option._
import cats.syntax.traverse._
import cats.syntax.validated._
import com.typesafe.config.ConfigFactory
import common.Checked
import common.validation.ErrorOr._
import common.validation.Validation._
import cwl.CommandLineTool._
import cwl.CwlType.CwlType
import cwl.CwlVersion._
import cwl.command.ParentName
import cwl.requirement.RequirementToAttributeMap
import eu.timepit.refined.W
import shapeless.syntax.singleton._
import shapeless.{:+:, CNil, Coproduct, Inl, Poly1, Witness}
import wom.callable.Callable.{InputDefinitionWithDefault, OutputDefinition, RequiredInputDefinition}
import wom.callable.{Callable, CallableTaskDefinition}
import wom.executable.Executable
import wom.expression.{ValueAsAnExpression, WomExpression}
import wom.types.{WomStringType, WomType}
import wom.values.{WomArray, WomEvaluatedCallInputs, WomFile, WomString, WomValue}
import wom.{CommandPart, RuntimeAttributes}

import scala.language.postfixOps
import scala.math.Ordering
import scala.util.Try

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
                                    permanentFailCodes: Option[Array[Int]]) {

  private [cwl] implicit val explicitWorkflowName = ParentName(id)
  private val inputNames = this.inputs.map(i => FullyQualifiedName(i.id).id).toSet

  // Circe can't create bidirectional links between workflow steps and runs (including `CommandLineTool`s) so this
  // ugly var is here to link back to a possible parent workflow step. This is needed to navigate upward for finding
  // requirements in the containment hierarchy. There isn't always a containing workflow step so this is an `Option`.
  private[cwl] var parentWorkflowStep: Option[WorkflowStep] = None

  /** Builds an `Executable` directly from a `CommandLineTool` CWL with no parent workflow. */
  def womExecutable(validator: RequirementsValidator, inputFile: Option[String] = None): Checked[Executable] = {
    val taskDefinition = buildTaskDefinition(validator, Vector.empty)
    CwlExecutableValidation.buildWomExecutable(taskDefinition, inputFile)
  }

  private def validateRequirementsAndHints(validator: RequirementsValidator): ErrorOr[List[Requirement]] = {
    import cats.instances.list._
    import cats.syntax.traverse._

    val allRequirements = requirements.toList.flatten ++ parentWorkflowStep.toList.flatMap(_.allRequirements)
    // All requirements must validate or this fails.
    val errorOrValidatedRequirements: ErrorOr[List[Requirement]] = allRequirements traverse validator

    errorOrValidatedRequirements map { validRequirements =>
      // Only Requirement hints, everything else is thrown out.
      // TODO CWL don't throw them out but pass them back to the caller to do with as the caller pleases.
      val hintRequirements = hints.toList.flatten.flatMap { _.select[Requirement] }
      val parentHintRequirements = parentWorkflowStep.toList.flatMap(_.allHints)

      // Throw out invalid Requirement hints.
      // TODO CWL pass invalid hints back to the caller to do with as the caller pleases.
      val validHints = (hintRequirements ++ parentHintRequirements).collect { case req if validator(req).isValid => req }
      validRequirements ++ validHints
    }
  }

  private def processRequirement(requirement: Requirement, expressionLib: ExpressionLib): Map[String, WomExpression] = {
    requirement.fold(RequirementToAttributeMap).apply(inputNames, expressionLib)
  }

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

        val defaultValue = inputParameter.default.map(_.fold(CommandInputParameter.DefaultToWomValuePoly).apply(womType))

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

  def buildTaskDefinition(validator: RequirementsValidator, parentExpressionLib: ExpressionLib): ErrorOr[CallableTaskDefinition] = {
    validateRequirementsAndHints(validator) map { requirementsAndHints: Seq[cwl.Requirement] =>
      val id = this.id

      val expressionLib: ExpressionLib =
        parentExpressionLib ++ inlineJavascriptRequirements(requirementsAndHints)

      // This is basically doing a `foldMap` but can't actually be a `foldMap` because:
      // - There is no monoid instance for `WomExpression`s.
      // - We want to fold from the right so the hints and requirements with the lowest precedence are processed first
      //   and later overridden if there are duplicate hints or requirements of the same type with higher precedence.
      val finalAttributesMap: Map[String, WomExpression] = (requirementsAndHints ++ DefaultDockerRequirement).foldRight(Map.empty[String, WomExpression])({
        case (requirement, attributesMap) => attributesMap ++ processRequirement(requirement, expressionLib)
      })

      val runtimeAttributes: RuntimeAttributes = RuntimeAttributes(finalAttributesMap)

      val meta: Map[String, String] = Map.empty
      val parameterMeta: Map[String, String] = Map.empty

      /*
    quoted from: http://www.commonwl.org/v1.0/CommandLineTool.html#CommandOutputBinding :

    For inputs and outputs, we only keep the variable name in the definition
     */

      val outputs: List[Callable.OutputDefinition] = this.outputs.map {
        case p @ CommandOutputParameter(cop_id, _, _, _, _, _, _, Some(tpe)) =>
          val womType = tpe.fold(MyriadOutputTypeToWomType)
          OutputDefinition(FullyQualifiedName(cop_id).id, womType, CommandOutputParameterExpression(p, womType, inputNames, expressionLib))
        case other => throw new NotImplementedError(s"Command output parameters such as $other are not yet supported")
      }.toList

    val inputDefinitions: List[_ <: Callable.InputDefinition] =
      this.inputs.map {
        case CommandInputParameter(inputId, _, _, _, _, _, _, Some(default), Some(tpe)) =>
          val inputType = tpe.fold(MyriadInputTypeToWomType)
          val inputName = FullyQualifiedName(inputId).id
          val defaultWomValue = default.fold(CommandInputParameter.DefaultToWomValuePoly).apply(inputType).toTry.get
          InputDefinitionWithDefault(inputName, inputType, ValueAsAnExpression(defaultWomValue))
        case CommandInputParameter(inputId, _, _, _, _, _, _, None, Some(tpe)) =>
          val inputType = tpe.fold(MyriadInputTypeToWomType)
          val inputName = FullyQualifiedName(inputId).id
          RequiredInputDefinition(inputName, inputType)
        case other => throw new NotImplementedError(s"command input parameters such as $other are not yet supported")
      }.toList

      def stringOrExpressionToString(soe: Option[StringOrExpression]): Option[String] = soe flatMap {
        case StringOrExpression.String(str) => Some(str)
        case StringOrExpression.Expression(_) => None // ... for now!
      }

      // The try will succeed if this is a task within a step. If it's a standalone file, the ID will be the file,
      // so the filename is the fallback.
      def taskName = Try(FullyQualifiedName(id).id).getOrElse(Paths.get(id).getFileName.toString)

      val adHocFileCreations: Set[WomExpression] = (for {
        requirements <- requirements.getOrElse(Array.empty[Requirement])
        initialWorkDirRequirement <- requirements.select[InitialWorkDirRequirement].toArray
        listing <- initialWorkDirRequirement.listings
      } yield InitialWorkDirFileGeneratorExpression(listing, expressionLib)).toSet[WomExpression]

      CallableTaskDefinition(
        taskName,
        buildCommandTemplate(expressionLib),
        runtimeAttributes,
        meta,
        parameterMeta,
        outputs,
        inputDefinitions,
        // TODO: This doesn't work in all cases and it feels clunky anyway - find a way to sort that out
        prefixSeparator = "#",
        commandPartSeparator = " ",
        stdoutRedirection = stringOrExpressionToString(stdout),
        stderrRedirection = stringOrExpressionToString(stderr),
        adHocFileCreation = adHocFileCreations
      )
    }
  }

  def asCwl = Coproduct[Cwl](this)

}

object CommandLineTool {
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

  def inlineJavascriptRequirements(allRequirementsAndHints: Seq[Requirement]): Vector[String] = {
    val inlineJavscriptRequirements: Seq[InlineJavascriptRequirement] = allRequirementsAndHints.toList.collect {
      case Inl(ijr:InlineJavascriptRequirement) => ijr
    }

    inlineJavscriptRequirements.flatMap(_.expressionLib.toList.flatten).toVector
  }

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
                                    secondaryFiles: Option[Array[StringOrExpression]] = None,
                                    format: Option[Expression :+: Array[String] :+: String :+: CNil] = None, //only valid when type: File
                                    streamable: Option[Boolean] = None, //only valid when type: File
                                    doc: Option[String :+: Array[String] :+: CNil] = None,
                                    inputBinding: Option[InputCommandLineBinding] = None,
                                    default: Option[CwlAny] = None,
                                    `type`: Option[MyriadInputType] = None)

  object CommandInputParameter {

    type DefaultToWomValueFunction = WomType => ErrorOr[WomValue]

    object DefaultToWomValuePoly extends Poly1 {
      implicit def caseFileOrDirectory: Case.Aux[FileOrDirectory, DefaultToWomValueFunction] = {
        at {
          _.fold(this)
        }
      }

      implicit def caseFileOrDirectoryArray: Case.Aux[Array[FileOrDirectory], DefaultToWomValueFunction] = {
        at {
          fileOrDirectoryArray =>
            womType =>
              fileOrDirectoryArray
                .toList
                .traverse[ErrorOr, WomValue](_.fold(this).apply(womType))
                .map(WomArray(_))
        }
      }

      implicit def caseFile: Case.Aux[File, DefaultToWomValueFunction] = {
        at {
          file =>
            womType =>
              file.asWomValue.flatMap(womType.coerceRawValue(_).toErrorOr)
        }
      }

      implicit def caseDirectory: Case.Aux[Directory, DefaultToWomValueFunction] = {
        at {
          directory =>
            womType =>
              directory.asWomValue.flatMap(womType.coerceRawValue(_).toErrorOr)
        }
      }

      implicit def caseJson: Case.Aux[io.circe.Json, DefaultToWomValueFunction] = {
        at {
          circeJson =>
            womType =>
              val stringJson = circeJson.noSpaces
              import spray.json._
              val sprayJson = stringJson.parseJson
              womType.coerceRawValue(sprayJson).toErrorOr
        }
      }
    }
  }

  case class CommandInputRecordSchema(
                                       `type`: W.`"record"`.T,
                                       fields: Option[Array[CommandInputRecordField]],
                                       label: Option[String])

  case class CommandInputRecordField(
                                      name: String,
                                      `type`: MyriadInputType,
                                      doc: Option[String],
                                      inputBinding: Option[InputCommandLineBinding],
                                      label: Option[String])

  case class CommandInputEnumSchema(
                                     symbols: Array[String],
                                     `type`: W.`"enum"`.T,
                                     label: Option[String],
                                     inputBinding: Option[InputCommandLineBinding])

  case class CommandInputArraySchema(
                                      items:
                                      CwlType :+:
                                        CommandInputRecordSchema :+:
                                        CommandInputEnumSchema :+:
                                        CommandInputArraySchema :+:
                                        String :+:
                                        Array[
                                          CwlType :+:
                                            CommandInputRecordSchema :+:
                                            CommandInputEnumSchema :+:
                                            CommandInputArraySchema :+:
                                            String :+:
                                            CNil] :+:
                                        CNil,
                                      `type`: W.`"array"`.T,
                                      label: Option[String],
                                      inputBinding: Option[InputCommandLineBinding])


  case class CommandOutputParameter(
                                     id: String,
                                     label: Option[String] = None,
                                     secondaryFiles: Option[SecondaryFiles] = None,
                                     format: Option[StringOrExpression] = None, //only valid when type: File
                                     streamable: Option[Boolean] = None, //only valid when type: File
                                     doc: Option[String :+: Array[String] :+: CNil] = None,
                                     outputBinding: Option[CommandOutputBinding] = None,
                                     `type`: Option[MyriadOutputType] = None)

  object CommandOutputParameter {

    def format(formatOption: Option[StringOrExpression],
               parameterContext: ParameterContext,
               expressionLib: ExpressionLib): ErrorOr[Option[String]] = {
      formatOption.traverse[ErrorOr, String] {
        format(_, parameterContext, expressionLib)
      }
    }

    def format(format: StringOrExpression, parameterContext: ParameterContext, expressionLib: ExpressionLib): ErrorOr[String] = {
      format.fold(CommandLineTool.CommandOutputParameter.FormatPoly).apply(parameterContext, expressionLib)
    }

    type FormatFunction = (ParameterContext, ExpressionLib) => ErrorOr[String]

    object FormatPoly extends Poly1 {
      implicit def caseStringOrExpression: Case.Aux[StringOrExpression, FormatFunction] = {
        at {
          _.fold(this)
        }
      }

      implicit def caseExpression: Case.Aux[Expression, FormatFunction] = {
        at {
          expression =>
            (parameterContext, expressionLib) =>
              val result: ErrorOr[WomValue] = ExpressionEvaluator.eval(expression, parameterContext, expressionLib)
              result flatMap {
                case womString: WomString => womString.value.valid
                case other => s"Not a valid file format: $other".invalidNel
              }
        }
      }

      implicit def caseString: Case.Aux[String, FormatFunction] = at { string => (_,_) => string.valid }
    }

    /**
      * Returns the list of secondary files for the primary file.
      */
    def secondaryFiles(primaryWomFile: WomFile,
                       secondaryFilesOption: Option[SecondaryFiles],
                       parameterContext: ParameterContext,
                       expressionLib: ExpressionLib): ErrorOr[List[WomFile]] = {
      secondaryFilesOption
        .map(secondaryFiles(primaryWomFile, _, parameterContext, expressionLib))
        .getOrElse(Nil.valid)
    }

    /**
      * Returns the list of secondary files for the primary file.
      */
    def secondaryFiles(primaryWomFile: WomFile,
                       secondaryFiles: SecondaryFiles,
                       parameterContext: ParameterContext,
                       expressionLib: ExpressionLib): ErrorOr[List[WomFile]] = {
      secondaryFiles
        .fold(CommandLineTool.CommandOutputParameter.SecondaryFilesPoly)
        .apply(primaryWomFile, parameterContext, expressionLib)
    }

    type SecondaryFilesFunction = (WomFile, ParameterContext, ExpressionLib) => ErrorOr[List[WomFile]]

    object SecondaryFilesPoly extends Poly1 {
      implicit def caseStringOrExpression: Case.Aux[StringOrExpression, SecondaryFilesFunction] = {
        at {
          _.fold(this)
        }
      }

      implicit def caseExpression: Case.Aux[Expression, SecondaryFilesFunction] = {
        at {
          expression =>
            (primaryWomFile, parameterContext, expressionLib) =>
              File.secondaryExpressionFiles(primaryWomFile, expression, parameterContext, expressionLib)
        }
      }

      implicit def caseString: Case.Aux[String, SecondaryFilesFunction] = {
        at {
          string =>
            (primaryWomFile, _, _) =>
              File.secondaryStringFile(primaryWomFile, string).map(List(_))
        }
      }

      implicit def caseArray: Case.Aux[Array[StringOrExpression], SecondaryFilesFunction] = {
        at {
          array =>
            (primaryWomFile, parameterContext, expressionLib) =>
              val functions: List[SecondaryFilesFunction] = array.toList.map(_.fold(this))
              functions.flatTraverse(_ (primaryWomFile, parameterContext, expressionLib))
        }
      }
    }

  }

  // Used to supply a default Docker image for platforms like PAPI that must have one even if the CWL document does
  // not specify a `DockerRequirement`.
  import net.ceedubs.ficus.Ficus._
  lazy val DefaultDockerImage = ConfigFactory.load().as[Option[String]]("cwl.default-docker-image")

  lazy val DefaultDockerRequirement = DefaultDockerImage map { image => Coproduct[Requirement](DockerRequirement(
    `class` = "DockerRequirement",
    dockerPull = Option(image),
    dockerLoad = None,
    dockerFile = None,
    dockerImport = None,
    dockerImageId = None,
    dockerOutputDirectory = None
  )) } toList

}
