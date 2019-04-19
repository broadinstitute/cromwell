package cwl

import java.nio.file.Paths

import common.Checked
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cwl.CwlVersion.CwlVersion
import cwl.Tool.inlineJavascriptRequirements
import cwl.command.ParentName
import cwl.requirement.RequirementToAttributeMap
import shapeless.Inl
import wom.callable.Callable.{OverridableInputDefinitionWithDefault, OptionalInputDefinition, OutputDefinition, RequiredInputDefinition}
import wom.callable.MetaValueElement.{MetaValueElementBoolean, MetaValueElementObject}
import wom.callable.{Callable, MetaValueElement, TaskDefinition}
import wom.executable.Executable
import wom.expression.{IoFunctionSet, ValueAsAnExpression, WomExpression}
import wom.types.WomOptionalType
import wom.values.{WomInteger, WomLong, WomString}
import wom.{RuntimeAttributes, RuntimeAttributesKeys}

import scala.util.Try

object Tool {
  def inlineJavascriptRequirements(allRequirementsAndHints: Seq[Requirement]): Vector[String] = {
    val inlineJavscriptRequirements: Seq[InlineJavascriptRequirement] = allRequirementsAndHints.toList.collect {
      case Inl(ijr:InlineJavascriptRequirement) => ijr
    }

    inlineJavscriptRequirements.flatMap(_.expressionLib.toList.flatten).toVector
  }
}

/**
  * Abstraction over a Tool (CommandLineTool or ExpressionTool)
  * Contains common logic to both, mostly to build the wom task definition.
  */
trait Tool {
  def inputs: Array[_ <: InputParameter]
  def outputs: Array[_ <: OutputParameter]
  def `class`: String
  def id: String
  def requirements: Option[Array[Requirement]]
  def hints: Option[Array[Hint]]
  def label: Option[String]
  def doc: Option[String]
  def cwlVersion: Option[CwlVersion]

  def asCwl: Cwl

  /** Builds an `Executable` directly from a `Tool` CWL with no parent workflow. */
  def womExecutable(validator: RequirementsValidator, inputFile: Option[String] = None, ioFunctions: IoFunctionSet, strictValidation: Boolean): Checked[Executable] = {
    val taskDefinition = buildTaskDefinition(validator, Vector.empty)
    CwlExecutableValidation.buildWomExecutable(taskDefinition, inputFile, ioFunctions, strictValidation)
  }

  protected def buildTaskDefinition(taskName: String,
                                    inputDefinitions: List[_ <: Callable.InputDefinition],
                                    outputDefinitions: List[Callable.OutputDefinition],
                                    runtimeAttributes: RuntimeAttributes,
                                    requirementsAndHints: List[cwl.Requirement],
                                    expressionLib: ExpressionLib): ErrorOr[TaskDefinition]

  private [cwl] implicit val explicitWorkflowName = ParentName(id)
  protected val inputNames: Set[String] = this.inputs.map(i => FullyQualifiedName(i.id).id).toSet

  // Circe can't create bidirectional links between workflow steps and runs (including `CommandLineTool`s) so this
  // ugly var is here to link back to a possible parent workflow step. This is needed to navigate upward for finding
  // requirements in the containment hierarchy. There isn't always a containing workflow step so this is an `Option`.
  private[cwl] var parentWorkflowStep: Option[WorkflowStep] = None


  private def validateRequirementsAndHints(validator: RequirementsValidator): ErrorOr[List[Requirement]] = {
    import cats.instances.list._
    import cats.syntax.traverse._

    val allRequirements = requirements.toList.flatten ++ parentWorkflowStep.toList.flatMap(_.allRequirements.list)
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

  protected def toolAttributes: Map[String, WomExpression] = Map.empty

  def buildTaskDefinition(validator: RequirementsValidator, parentExpressionLib: ExpressionLib): Checked[TaskDefinition] = {
    def build(requirementsAndHints: Seq[cwl.Requirement]) = {
      val id = this.id

      val expressionLib: ExpressionLib = parentExpressionLib ++ inlineJavascriptRequirements(requirementsAndHints)

      // This is basically doing a `foldMap` but can't actually be a `foldMap` because:
      // - There is no monoid instance for `WomExpression`s.
      // - We want to fold from the right so the hints and requirements with the lowest precedence are processed first
      //   and later overridden if there are duplicate hints or requirements of the same type with higher precedence.
      val attributesMap: Map[String, WomExpression] = requirementsAndHints.foldRight(Map.empty[String, WomExpression])({
        case (requirement, acc) => acc ++ processRequirement(requirement, expressionLib)
      })

      val runtimeAttributes: RuntimeAttributes = RuntimeAttributes(attributesMap ++ toolAttributes)

      val schemaDefRequirement: SchemaDefRequirement = requirementsAndHints.flatMap{
        _.select[SchemaDefRequirement].toList
      }.headOption.getOrElse(SchemaDefRequirement())

      lazy val localizationOptionalMetaObject: MetaValueElement = MetaValueElementObject(
        Map(
          "localization_optional" -> MetaValueElementBoolean(true)
        )
      )

      // If input dir min == 0 that's a signal not to localize the input files
      // https://github.com/dnanexus/dx-cwl/blob/fca163d825beb62f8a3004f2e0a6742805e6218c/dx-cwl#L426
      val localizationOptional = runtimeAttributes.attributes.get(RuntimeAttributesKeys.DnaNexusInputDirMinKey) match {
        case Some(ValueAsAnExpression(WomInteger(0))) => Option(localizationOptionalMetaObject)
        case Some(ValueAsAnExpression(WomLong(0L))) => Option(localizationOptionalMetaObject)
        case Some(ValueAsAnExpression(WomString("0"))) => Option(localizationOptionalMetaObject)
        case _ => None
      }

      val inputDefinitions: List[_ <: Callable.InputDefinition] =
        this.inputs.map {
          case input @ InputParameter.IdDefaultAndType(inputId, default, tpe) =>
            val inputType = tpe.fold(MyriadInputTypeToWomType).apply(schemaDefRequirement)
            val inputName = FullyQualifiedName(inputId).id
            val defaultWomValue = default.fold(InputParameter.DefaultToWomValuePoly).apply(inputType).toTry.get
            OverridableInputDefinitionWithDefault(
              inputName,
              inputType,
              ValueAsAnExpression(defaultWomValue),
              InputParameter.inputValueMapper(input, tpe, expressionLib, asCwl.schemaOption),
              localizationOptional
            )
          case input @ InputParameter.IdAndType(inputId, tpe) =>
            val inputType = tpe.fold(MyriadInputTypeToWomType).apply(schemaDefRequirement)
            val inputName = FullyQualifiedName(inputId).id
            inputType match {
              case optional: WomOptionalType =>
                OptionalInputDefinition(
                  inputName,
                  optional,
                  InputParameter.inputValueMapper(input, tpe, expressionLib, asCwl.schemaOption),
                  localizationOptional
                )
              case _ =>
                RequiredInputDefinition(
                  inputName,
                  inputType,
                  InputParameter.inputValueMapper(input, tpe, expressionLib, asCwl.schemaOption),
                  localizationOptional
                )
            }
          case other => throw new UnsupportedOperationException(s"command input parameters such as $other are not yet supported")
        }.toList

      val outputDefinitions: List[Callable.OutputDefinition] = this.outputs.map {
        case p @ OutputParameter.IdAndType(cop_id, tpe) =>
          val womType = tpe.fold(MyriadOutputTypeToWomType).apply(schemaDefRequirement)
          OutputDefinition(FullyQualifiedName(cop_id).id, womType, OutputParameterExpression(p, womType, inputNames, expressionLib, schemaDefRequirement))
        case other => throw new UnsupportedOperationException(s"Command output parameters such as $other are not yet supported")
      }.toList

      // The try will succeed if this is a task within a step. If it's a standalone file, the ID will be the file,
      // so the filename is the fallback.
      def taskName = Try(FullyQualifiedName(id).id).getOrElse(Paths.get(id).getFileName.toString)

      buildTaskDefinition(taskName, inputDefinitions, outputDefinitions, runtimeAttributes, requirementsAndHints.toList, expressionLib).toEither
    }

    validateRequirementsAndHints(validator).toEither match {
      case Right(requirementsAndHints: Seq[cwl.Requirement]) => build(requirementsAndHints)
      case Left(errors) => Left(errors)
    }
  }
}
