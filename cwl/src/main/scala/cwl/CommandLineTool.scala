package cwl

import java.nio.file.Paths

import common.Checked
import common.validation.ErrorOr.ErrorOr
import cwl.CommandLineTool._
import cwl.CwlType.CwlType
import cwl.CwlVersion._
import cwl.command.ParentName
import cwl.requirement.ResourceRequirementToWomExpression
import eu.timepit.refined.W
import shapeless.syntax.singleton._
import shapeless.{:+:, CNil, Coproduct, Witness}
import wom.callable.Callable.{InputDefinitionWithDefault, OutputDefinition, RequiredInputDefinition}
import wom.callable.{Callable, CallableTaskDefinition}
import wom.executable.Executable
import wom.expression.{InputLookupExpression, ValueAsAnExpression, WomExpression}
import wom.types.WomType
import wom.values.WomString
import wom.{CommandPart, RuntimeAttributes}

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

  /** Builds an `Executable` directly from a `CommandLineTool` CWL with no parent workflow. */
  def womExecutable(validator: RequirementsValidator, inputFile: Option[String] = None): Checked[Executable] = {
    val taskDefinition = buildTaskDefinition(parentWorkflow = None, validator)
    CwlExecutableValidation.buildWomExecutable(taskDefinition, inputFile)
  }

  private def validateRequirementsAndHints(parentWorkflow: Option[Workflow], validator: RequirementsValidator): ErrorOr[List[Requirement]] = {
    import cats.instances.list._
    import cats.syntax.traverse._

    val allRequirements = requirements.toList.flatten ++ parentWorkflow.toList.flatMap(_.allRequirements)
    // All requirements must validate or this fails.
    val errorOrValidatedRequirements: ErrorOr[List[Requirement]] = allRequirements traverse validator

    errorOrValidatedRequirements map { validRequirements =>
      // Only Requirement hints, everything else is thrown out.
      // TODO CWL don't throw them out but pass them back to the caller to do with as the caller pleases.
      val hintRequirements = hints.toList.flatten.flatMap { _.select[Requirement] }
      val parentHintRequirements = parentWorkflow.toList.flatMap(_.allHints)

      // Throw out invalid Requirement hints.
      // TODO CWL pass invalid hints back to the caller to do with as the caller pleases.
      val validHints = (hintRequirements ++ parentHintRequirements).collect { case req if validator(req).isValid => req }
      validRequirements ++ validHints
    }
  }

  def buildTaskDefinition(parentWorkflow: Option[Workflow], validator: RequirementsValidator): ErrorOr[CallableTaskDefinition] = {
    import wom.RuntimeAttributesKeys._
    
    def processRequirement(requirement: Requirement): Map[String, WomExpression] = {
      def processDockerRequirement(requirement: DockerRequirement): Map[String, WomExpression] = {
        requirement.dockerPull.orElse(requirement.dockerImageId).map({ pull =>
          DockerKey -> ValueAsAnExpression(WomString(pull))
        }).toMap
      }

      def processResourceRequirement(requirement: ResourceRequirement): Map[String, WomExpression] = {
        def toExpression(resourceRequirement: ResourceRequirementType) = resourceRequirement.fold(ResourceRequirementToWomExpression).apply(inputNames)

        List(
          // Map cpuMin to both cpuMin and cpu keys
          requirement.effectiveCoreMin.toList.map(toExpression).flatMap(min => List(CPUMinKey -> min, CPUKey -> min)),
          requirement.effectiveCoreMax.toList.map(toExpression).map(CPUMaxKey -> _),
          // Map ramMin to both memoryMin and memory keys
          requirement.effectiveRamMin.toList.map(toExpression).flatMap(min => List(MemoryMinKey -> min, MemoryKey -> min)),
          requirement.effectiveRamMax.toList.map(toExpression).map(MemoryMaxKey -> _),
          requirement.effectiveTmpdirMin.toList.map(toExpression).map(TmpDirMinKey -> _),
          requirement.effectiveTmpdirMax.toList.map(toExpression).map(TmpDirMaxKey -> _),
          requirement.effectiveOutdirMin.toList.map(toExpression).map(OutDirMinKey -> _),
          requirement.effectiveOutdirMax.toList.map(toExpression).map(OutDirMaxKey -> _)
        ).flatten.toMap
      }

      requirement.select[DockerRequirement].map(processDockerRequirement).orElse(
        requirement.select[ResourceRequirement].map(processResourceRequirement)
      ).getOrElse(Map.empty)
    }
    
    validateRequirementsAndHints(parentWorkflow, validator) map { requirementsAndHints =>
      val id = this.id

      val commandTemplate: Seq[CommandPart] = baseCommand.toSeq.flatMap(_.fold(BaseCommandToCommandParts)) ++
        arguments.toSeq.flatMap(_.map(_.fold(ArgumentToCommandPart))) ++
        CommandLineTool.orderedForCommandLine(inputs).map(InputParameterCommandPart.apply)

      val finalAttributesMap = requirementsAndHints.foldRight(Map.empty[String, WomExpression])({
        case (requirement, attributesMap) => attributesMap ++ processRequirement(requirement)
      })

      val runtimeAttributes: RuntimeAttributes = RuntimeAttributes(finalAttributesMap)

      val meta: Map[String, String] = Map.empty
      val parameterMeta: Map[String, String] = Map.empty

      /*
    quoted from: http://www.commonwl.org/v1.0/CommandLineTool.html#CommandOutputBinding :

    For inputs and outputs, we only keep the variable name in the definition

    The output parameter value is generated by applying these operations in the following order:

      glob
      loadContents
      outputEval
      secondaryFiles
     */

    val outputs: List[Callable.OutputDefinition] = this.outputs.map {
      case CommandOutputParameter(cop_id, _, _, _, _, _, Some(outputBinding), Some(tpe)) =>
        val womType = tpe.fold(MyriadOutputTypeToWomType)
        OutputDefinition(FullyQualifiedName(cop_id).id, womType, CommandOutputExpression(outputBinding, womType, inputNames))

      //This catches states where the output binding is not declared but the type is
      case CommandOutputParameter(id, _, _, _, _, _, _, Some(tpe)) =>
        val womType: WomType = tpe.fold(MyriadOutputTypeToWomType)
        OutputDefinition(FullyQualifiedName(id).id, womType, InputLookupExpression(womType, id))

      case other => throw new NotImplementedError(s"Command output parameters such as $other are not yet supported")
    }.toList

    val inputDefinitions: List[_ <: Callable.InputDefinition] =
      this.inputs.map {
        case CommandInputParameter(id, _, _, _, _, _, _, Some(default), Some(tpe)) =>
          val inputType = tpe.fold(MyriadInputTypeToWomType)
          val inputName = FullyQualifiedName(id).id
          InputDefinitionWithDefault(inputName, inputType, ValueAsAnExpression(inputType.coerceRawValue(default.toString).get))
        case CommandInputParameter(id, _, _, _, _, _, _, None, Some(tpe)) =>
          val inputType = tpe.fold(MyriadInputTypeToWomType)
          val inputName = FullyQualifiedName(id).id
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
      } yield InitialWorkDirFileGeneratorExpression(listing)).toSet[WomExpression]

      CallableTaskDefinition(
        taskName,
        commandTemplate,
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

  /**
    * Sort according to position. If position does not exist, use 0 per spec:
    * http://www.commonwl.org/v1.0/CommandLineTool.html#CommandLineBinding
    *
    * If an input binding is not specified, ignore the input parameter.
    */
  protected[cwl] def orderedForCommandLine(inputs: Array[CommandInputParameter]): Seq[CommandInputParameter] = {
    inputs.
      filter(_.inputBinding.isDefined).
      sortBy(_.inputBinding.flatMap(_.position).getOrElse(0)).
      toSeq
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

  type BaseCommand = String :+: Array[String] :+: CNil

  type Argument = Expression :+: CommandLineBinding :+: String :+: CNil

  case class CommandInputParameter(
                                    id: String,
                                    label: Option[String] = None,
                                    secondaryFiles: Option[Array[Expression :+: String :+: CNil]] = None,
                                    format: Option[Expression :+: Array[String] :+: String :+: CNil] = None, //only valid when type: File
                                    streamable: Option[Boolean] = None, //only valid when type: File
                                    doc: Option[String :+: Array[String] :+: CNil] = None,
                                    inputBinding: Option[CommandLineBinding] = None,
                                    default: Option[CwlAny] = None,
                                    `type`: Option[MyriadInputType] = None)

  case class CommandInputRecordSchema(
                                       `type`: W.`"record"`.T,
                                       fields: Option[Array[CommandInputRecordField]],
                                       label: Option[String])

  case class CommandInputRecordField(
                                      name: String,
                                      `type`: MyriadInputType,
                                      doc: Option[String],
                                      inputBinding: Option[CommandLineBinding],
                                      label: Option[String])

  case class CommandInputEnumSchema(
                                     symbols: Array[String],
                                     `type`: W.`"enum"`.T,
                                     label: Option[String],
                                     inputBinding: Option[CommandLineBinding])

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
                                      inputBinding: Option[CommandLineBinding])


  case class CommandOutputParameter(
                                     id: String,
                                     label: Option[String] = None,
                                     secondaryFiles: Option[Expression :+: String :+: Array[Expression :+: String :+: CNil] :+: CNil] = None,
                                     format: Option[Expression :+: Array[String] :+: String :+: CNil] = None, //only valid when type: File
                                     streamable: Option[Boolean] = None, //only valid when type: File
                                     doc: Option[String :+: Array[String] :+: CNil] = None,
                                     outputBinding: Option[CommandOutputBinding] = None,
                                     `type`: Option[MyriadOutputType] = None)


}
