package wdl.transforms.v1_2.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.collections.EnhancedCollections._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.Checked._
import common.validation.ErrorOr._
import wdl.model.draft3.elements.CommandPartElement.{PlaceholderCommandPartElement, StringCommandPartElement}
import wdl.model.draft3.elements.ExpressionElement.{IdentifierLookup, KvPair, ObjectLiteral}
import wdl.model.draft3.elements._
import wdl.model.v1_2.elements.V1_2ExpressionElement.EnvironmentVariableValue
import wdl.model.v1_2.elements.{HintsSectionElement, RequirementsSectionElement, WorkflowHintsSectionElement}
import wdl.transforms.base.ast2wdlom._
import wom.SourceFileLocation
import wom.callable.MetaKvPair
import wom.callable.MetaValueElement.MetaValueElementBoolean

/**
 * AST conversions for WDL 1.2 language features. Each converter handles the new syntax and delegates everything else
 * to the shared WDL 1.0+ converters in `wdl.transforms.base.ast2wdlom`.
 */
object AstToV1_2Elements {

  private val hintsObjectLiteralNames = Set("HintsObjectLiteral", "InputHintsObjectLiteral", "OutputHintsObjectLiteral")

  // Task sections

  def astToTaskSectionElement(baseConverter: CheckedAtoB[GenericAst, TaskSectionElement])(implicit
    astNodeToKvPair: CheckedAtoB[GenericAstNode, KvPair],
    astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement]
  ): CheckedAtoB[GenericAst, TaskSectionElement] = CheckedAtoB.fromCheck { a: GenericAst =>
    a.getName match {
      case "Requirements" => a.getAttributeAsVector[KvPair]("map").map(RequirementsSectionElement)
      case "Hints" => a.getAttributeAsVectorF("map")(astNodeToHintKvPair).map(HintsSectionElement)
      case _ => baseConverter.run(a)
    }
  }

  // A task hint value may be any expression, or one of the hints-scoped `hints`, `input` or `output` literals
  private def astNodeToHintKvPair(node: GenericAstNode)(implicit
    astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement]
  ): Checked[KvPair] = node match {
    case a: GenericAst if a.getName == "RuntimeAttribute" =>
      for {
        key <- a.getAttributeAs[String]("key")
        value <- astNodeToHintValue(a.getAttribute("value"))
      } yield KvPair(key, value)
    case other => s"Unexpected element in hints section: $other".invalidNelCheck
  }

  private def astNodeToHintValue(node: GenericAstNode)(implicit
    astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement]
  ): Checked[ExpressionElement] = node match {
    case a: GenericAst if hintsObjectLiteralNames.contains(a.getName) =>
      a.getAttributeAsVectorF("map")(astNodeToHintsObjectMember) flatMap { members =>
        val duplicates = members.groupBy(_._1).collect { case (key, values) if values.size > 1 => key }
        if (duplicates.isEmpty) ObjectLiteral(members.toMap).validNelCheck
        else s"Duplicate keys in hints object${a.lineAndColumnString}: ${duplicates.mkString(", ")}".invalidNelCheck
      }
    case other => astNodeToExpressionElement.run(other)
  }

  // Keys of `input` and `output` hints literals may use dotted notation to refer to struct members, e.g. `person.cv`
  private def astNodeToHintsObjectMember(node: GenericAstNode)(implicit
    astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement]
  ): Checked[(String, ExpressionElement)] = node match {
    case a: GenericAst if a.getName == "HintsObjectKV" =>
      for {
        key <- a.getAttributeAs[String]("key")
        members <- a.getAttributeAsVector[String]("members", optional = true)
        value <- astNodeToHintValue(a.getAttribute("value"))
      } yield (key +: members).mkString(".") -> value
    case other => s"Unexpected element in hints object: $other".invalidNelCheck
  }

  // Tasks

  def astToTaskDefinitionElement(implicit
    astNodeToTaskSectionElement: CheckedAtoB[GenericAstNode, TaskSectionElement]
  ): CheckedAtoB[GenericAst, TaskDefinitionElement] = CheckedAtoB.fromErrorOr { a: GenericAst =>
    s"read task definition${a.lineAndColumnString}"
  } { a =>
    val nameValidation: ErrorOr[String] = astNodeToString(a.getAttribute("name")).toValidated
    val sectionsValidation: ErrorOr[Vector[TaskSectionElement]] =
      a.getAttributeAsVector[TaskSectionElement]("sections").toValidated
    val envNamesValidation: ErrorOr[Vector[String]] = envDeclarationNames(a).toValidated
    val sourceLocation: Option[SourceFileLocation] = a.getSourceLine.map(SourceFileLocation(_))

    (nameValidation, sectionsValidation, envNamesValidation) flatMapN { (name, sections, envNames) =>
      for {
        runtimeSections <- convertRequirementsAndHints(sections)
        task <- AstToTaskDefinitionElement.combineElements(name, runtimeSections, sourceLocation)
        withTaskVariable <- TaskRuntimeVariable.addIfReferenced(task)
      } yield exportEnvironmentVariables(withTaskVariable, envNames)
    }
  }

  /**
   * The `runtime` section is deprecated in WDL 1.2 and is mutually exclusive with `requirements` and `hints`. Either
   * way, the result is a single runtime attributes section so that the shared WOM conversion can be reused.
   */
  private def convertRequirementsAndHints(sections: Vector[TaskSectionElement]): ErrorOr[Vector[TaskSectionElement]] = {
    val runtimeSections = sections.filterByType[RuntimeAttributesSectionElement]
    val requirementsSections = sections.filterByType[RequirementsSectionElement]
    val hintsSections = sections.filterByType[HintsSectionElement]
    val otherSections = sections filter {
      case _: RuntimeAttributesSectionElement | _: RequirementsSectionElement | _: HintsSectionElement => false
      case _ => true
    }

    if (runtimeSections.nonEmpty && (requirementsSections.nonEmpty || hintsSections.nonEmpty)) {
      "A task cannot have a runtime section together with a requirements or hints section. The runtime section is deprecated in WDL 1.2; use requirements and hints instead.".invalidNel
    } else if (requirementsSections.size > 1) {
      s"Task must have no more than one requirements section, found ${requirementsSections.size}.".invalidNel
    } else if (hintsSections.size > 1) {
      s"Task must have no more than one hints section, found ${hintsSections.size}.".invalidNel
    } else if (runtimeSections.nonEmpty) {
      // TODO: The runtime section is deprecated in WDL 1.2. Suggestion: surface a deprecation warning in womtool
      //  validate and in the workflow's metadata once there is a mechanism for language factories to report warnings.
      runtimeSections.traverse(RequirementsAndHints.normalizeRuntimeSection).map(otherSections ++ _)
    } else {
      RequirementsAndHints.toRuntimeAttributesSection(requirementsSections.headOption, hintsSections.headOption) map {
        otherSections ++ _
      }
    }
  }

  // Names of input and private declarations annotated with the `env` modifier, in declaration order
  private def envDeclarationNames(task: GenericAst): Checked[Vector[String]] = {
    def envName(declaration: GenericAstNode): Checked[Option[String]] = declaration match {
      case d: GenericAst if Option(d.getAttribute("env")).isDefined =>
        d.getAttributeAs[String]("name").map(Option.apply)
      case _ => None.validNelCheck
    }

    task.getAttributeAsVectorF("sections")(_.validNelCheck) flatMap { sections =>
      sections.flatTraverse {
        case s: GenericAst if s.getName == "Inputs" =>
          s.getAttributeAsVectorF("inputs", optional = true)(envName).map(_.flatten)
        case s: GenericAst if s.getName == "Declaration" => envName(s).map(_.toVector)
        case _ => Vector.empty[String].validNelCheck
      }
    }
  }

  /**
   * Makes declarations annotated with `env` available as environment variables by exporting them at the start of the
   * command. The value is serialized and single-quoted, so it is never interpreted by the shell.
   *
   * https://github.com/openwdl/wdl/blob/v1.2.1/SPEC.md#environment-variables
   *
   * TODO: The spec leaves the environment variable length and count limits to the engine. Suggestion: document the
   *  limits of each backend (e.g. the GCP Batch and AWS Batch script size limits) in docs/LanguageSupport.md.
   */
  private def exportEnvironmentVariables(task: TaskDefinitionElement, envNames: Vector[String]): TaskDefinitionElement =
    if (envNames.isEmpty) task
    else {
      val exports = envNames map { name =>
        CommandSectionLine(
          Vector(
            StringCommandPartElement(s"export $name="),
            PlaceholderCommandPartElement(EnvironmentVariableValue(IdentifierLookup(name)),
                                          PlaceholderAttributeSet.empty
            )
          )
        )
      }
      task.copy(commandSection = CommandSectionElement(exports ++ task.commandSection.parts))
    }

  // Workflows

  def astToWorkflowBodyElement(baseConverter: CheckedAtoB[GenericAst, WorkflowBodyElement])(implicit
    astNodeToMetaKvPair: CheckedAtoB[GenericAstNode, MetaKvPair]
  ): CheckedAtoB[GenericAst, WorkflowBodyElement] = CheckedAtoB.fromCheck { a: GenericAst =>
    a.getName match {
      case "WorkflowHints" => a.getAttributeAsVector[MetaKvPair]("map").map(WorkflowHintsSectionElement)
      case _ => baseConverter.run(a)
    }
  }

  def astToWorkflowDefinitionElement(baseConverter: CheckedAtoB[GenericAst, WorkflowDefinitionElement])(implicit
    astNodeToWorkflowBodyElement: CheckedAtoB[GenericAstNode, WorkflowBodyElement]
  ): CheckedAtoB[GenericAst, WorkflowDefinitionElement] = CheckedAtoB.fromCheck { a: GenericAst =>
    for {
      _ <- validateNoEnvModifier(a)
      bodyElements <- a.getAttributeAsVector[WorkflowBodyElement]("body")
      hints <- bodyElements.filterByType[WorkflowHintsSectionElement] match {
        case Vector() => None.validNelCheck
        case Vector(single) => Option(single).validNelCheck
        case multiple => s"Workflow must have no more than one hints section, found ${multiple.size}.".invalidNelCheck
      }
      workflow <- baseConverter.run(a)
    } yield applyWorkflowHints(workflow, hints)
  }

  // The `env` modifier is only meaningful for task declarations
  private def validateNoEnvModifier(workflow: GenericAst): Checked[Unit] = {
    def containsEnvModifier(node: GenericAstNode): Boolean = node match {
      case a: GenericAst => a.getName == "EnvModifier" || a.getAttributes.values.exists(containsEnvModifier)
      case l: GenericAstList => l.astNodeList.exists(containsEnvModifier)
      case _ => false
    }
    if (containsEnvModifier(workflow))
      s"The env modifier is only allowed on task declarations, not in workflow${workflow.lineAndColumnString}".invalidNelCheck
    else ().validNelCheck
  }

  /**
   * In WDL 1.2, `allow_nested_inputs` (alias `allowNestedInputs`) moved from the workflow `meta` section to the
   * workflow `hints` section. It is copied into the meta section, where the shared WOM conversion reads it from.
   *
   * TODO: WDL 1.2 deprecates using `allow_nested_inputs` to leave *required* call inputs unsatisfied: only call inputs
   *  that have a default value may be set at runtime. Cromwell still allows required nested inputs, as in WDL 1.1.
   *  Suggestion: decide whether to enforce the stricter rule (a breaking change for workflows upgraded from 1.1) or to
   *  warn, and implement it in `WorkflowDefinitionElementToWomWorkflowDefinition` behind a version-specific flag.
   */
  private def applyWorkflowHints(workflow: WorkflowDefinitionElement,
                                 hints: Option[WorkflowHintsSectionElement]
  ): WorkflowDefinitionElement = {
    val allowNestedInputs = hints.toVector.flatMap(_.hints) collectFirst {
      case MetaKvPair("allow_nested_inputs" | "allowNestedInputs", value: MetaValueElementBoolean) => value
    }
    allowNestedInputs match {
      case Some(value) =>
        val existingMeta = workflow.metaSection.map(_.meta).getOrElse(Map.empty)
        workflow.copy(metaSection = Option(MetaSectionElement(existingMeta + ("allowNestedInputs" -> value))))
      case None => workflow
    }
  }

  // Structs

  /**
   * WDL 1.2 structs may contain `meta` and `parameter_meta` sections. They only document the struct, so they are
   * validated and then dropped.
   */
  def astToStructElement(implicit
    astNodeToTypeElement: CheckedAtoB[GenericAstNode, TypeElement],
    astNodeToMetaSectionElement: CheckedAtoB[GenericAstNode, MetaSectionElement],
    astNodeToParameterMetaSectionElement: CheckedAtoB[GenericAstNode, ParameterMetaSectionElement]
  ): CheckedAtoB[GenericAst, StructElement] = CheckedAtoB.fromErrorOr("convert AST to struct definition") {
    a: GenericAst =>
      val nameValidation: ErrorOr[String] = a.getAttributeAs[String]("name").toValidated
      val elementsValidation: ErrorOr[Vector[GenericAstNode]] =
        a.getAttributeAsVectorF("entries")(_.validNelCheck).toValidated

      (nameValidation, elementsValidation) flatMapN { (name, elements) =>
        val entriesValidation: ErrorOr[Vector[StructEntryElement]] = elements.collect {
          case e: GenericAst if e.getName == "StructEntry" =>
            ((e.getAttributeAs[String]("name").toValidated: ErrorOr[String]),
             e.getAttributeAs[TypeElement]("type").toValidated: ErrorOr[TypeElement]
            ).mapN(StructEntryElement.apply)
        }.sequence
        val metaValidation: ErrorOr[Vector[MetaSectionElement]] =
          elements
            .collect { case e: GenericAst if e.getName == "Meta" => astNodeToMetaSectionElement.run(e) }
            .traverse(_.toValidated)
        val parameterMetaValidation: ErrorOr[Vector[ParameterMetaSectionElement]] =
          elements
            .collect {
              case e: GenericAst if e.getName == "ParameterMeta" => astNodeToParameterMetaSectionElement.run(e)
            }
            .traverse(_.toValidated)

        (entriesValidation, metaValidation, parameterMetaValidation) flatMapN { (entries, meta, parameterMeta) =>
          val memberNames = entries.map(_.identifier).toSet
          val unknownParameterMetaKeys = parameterMeta.flatMap(_.metaAttributes.keys).filterNot(memberNames.contains)
          if (meta.size > 1) s"Struct $name must have no more than one meta section, found ${meta.size}.".invalidNel
          else if (parameterMeta.size > 1)
            s"Struct $name must have no more than one parameter_meta section, found ${parameterMeta.size}.".invalidNel
          else if (unknownParameterMetaKeys.nonEmpty)
            s"Struct $name parameter_meta refers to members that do not exist: ${unknownParameterMetaKeys.mkString(", ")}".invalidNel
          else StructElement(name, entries).validNel
        }
      }
  }
}
