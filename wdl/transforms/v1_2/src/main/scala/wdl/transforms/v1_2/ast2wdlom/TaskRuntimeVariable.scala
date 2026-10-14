package wdl.transforms.v1_2.ast2wdlom

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.CommandPartElement.PlaceholderCommandPartElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.elements._
import wdl.model.draft3.graph.{
  UnlinkedCallOutputOrIdentifierAndMemberAccessHook,
  UnlinkedConsumedValueHook,
  UnlinkedIdentifierHook
}
import wdl.transforms.v1_2.linking.expression.consumed.expressionElementUnlinkedValueConsumer
import wom.RuntimeAttributesKeys
import wom.callable.MetaValueElement
import wom.callable.MetaValueElement._
import wom.values.{WomBoolean, WomFloat, WomInteger}

/**
 * The hidden `task` variable of WDL 1.2, which gives the `command` and `output` sections access to the task's name,
 * requirements and metadata at runtime.
 *
 * https://github.com/openwdl/wdl/blob/v1.2.1/SPEC.md#runtime-access-to-requirements-hints-and-metadata
 *
 * When a task refers to `task`, a private declaration of that name is added to the task. Its value is computed from
 * the task definition: per the spec, "if the runtime engine is not able to provide the actual value of a requirement,
 * then it must provide the requested value instead, or the default value if no specific value was requested".
 */
object TaskRuntimeVariable {

  val Name = "task"

  /**
   * Members that Cromwell cannot provide yet. Referring to them is a validation error rather than a wrong value.
   *
   * TODO: Implement these members. Suggestions:
   *  - `attempt`: needs the job's attempt number, which is only known to the engine. Inject it as a hidden input when
   *    the engine evaluates the task's declarations (`CallPreparation`/`JobPreparationActor`).
   *  - `return_code` (output section only): needs the command's return code. The engine could provide it from the job's
   *    `rc` file when evaluating outputs.
   *  - `memory`: is an Int number of bytes, but `WomInteger` is 32-bit, so the default of 2 GiB does not fit. WDL 1.2
   *    requires 64-bit Ints; changing `WomInteger` to a Long affects every WDL version.
   *  - `gpu`, `fpga`, `disks`: need the resources actually allocated by the backend.
   */
  val UnsupportedMembers: Set[String] = Set("attempt", "return_code", "memory", "gpu", "fpga", "disks")

  def addIfReferenced(task: TaskDefinitionElement): ErrorOr[TaskDefinitionElement] = {
    val commandAndOutputReferences = taskReferences(commandExpressions(task) ++ outputExpressions(task))
    val otherReferences = taskReferences(otherExpressions(task))

    if (otherReferences.nonEmpty) {
      s"The '$Name' variable is only available in the command and output sections of task '${task.name}'".invalidNel
    } else if (commandAndOutputReferences.isEmpty) {
      task.validNel
    } else {
      val unsupported = commandAndOutputReferences.collect {
        case UnlinkedCallOutputOrIdentifierAndMemberAccessHook(_, member) if UnsupportedMembers.contains(member) =>
          s"$Name.$member"
      }
      if (unsupported.nonEmpty)
        s"Task '${task.name}' uses ${unsupported.toList.sorted.mkString(", ")}, which Cromwell does not support yet".invalidNel
      else
        task.copy(declarations = task.declarations :+ declaration(task)).validNel
    }
  }

  private def declaration(task: TaskDefinitionElement): IntermediateValueDeclarationElement = {
    val runtime: Map[String, ExpressionElement] =
      task.runtimeSection.toVector.flatMap(_.runtimeAttributes).map(kv => kv.key -> kv.value).toMap

    val container = runtime
      .get(RuntimeAttributesKeys.ContainerKey)
      .orElse(runtime.get(RuntimeAttributesKeys.DockerKey))
      .getOrElse(NoneLiteralElement)

    val members: Map[String, ExpressionElement] = Map(
      "name" -> StringLiteral(task.name),
      // TODO: The spec suggests that the ID includes the call alias and the scatter index, which are only known to the
      //  engine. Suggestion: provide it together with `attempt` (see `UnsupportedMembers`).
      "id" -> StringLiteral(task.name),
      "container" -> container,
      "cpu" -> runtime.getOrElse(RuntimeAttributesKeys.CpuKey, PrimitiveLiteralExpressionElement(WomFloat(1.0))),
      // `None` means that the engine cannot determine whether the runtime of the task is limited
      "end_time" -> NoneLiteralElement,
      "meta" -> ObjectLiteral(task.metaSection.map(_.meta).getOrElse(Map.empty).map { case (k, v) =>
        k -> metaValueToExpression(v)
      }),
      "parameter_meta" -> ObjectLiteral(
        task.parameterMetaSection.map(_.metaAttributes).getOrElse(Map.empty).map { case (k, v) =>
          k -> metaValueToExpression(v)
        }
      ),
      "ext" -> ObjectLiteral(Map.empty)
    )

    IntermediateValueDeclarationElement(ObjectTypeElement, Name, ObjectLiteral(members))
  }

  private def metaValueToExpression(value: MetaValueElement): ExpressionElement = value match {
    case MetaValueElementNull => NoneLiteralElement
    case MetaValueElementBoolean(b) => PrimitiveLiteralExpressionElement(WomBoolean(b))
    case MetaValueElementFloat(f) => PrimitiveLiteralExpressionElement(WomFloat(f))
    case MetaValueElementInteger(i) => PrimitiveLiteralExpressionElement(WomInteger(i))
    case MetaValueElementString(s) => StringLiteral(s)
    case MetaValueElementObject(values) => ObjectLiteral(values.map { case (k, v) => k -> metaValueToExpression(v) })
    case MetaValueElementArray(values) => ArrayLiteral(values.map(metaValueToExpression))
  }

  private def taskReferences(expressions: Seq[ExpressionElement]): Set[UnlinkedConsumedValueHook] =
    expressions.toSet[ExpressionElement].flatMap { e =>
      expressionElementUnlinkedValueConsumer.expressionConsumedValueHooks(e)(expressionElementUnlinkedValueConsumer)
    } filter {
      case UnlinkedIdentifierHook(Name) => true
      case UnlinkedCallOutputOrIdentifierAndMemberAccessHook(Name, _) => true
      case _ => false
    }

  private def commandExpressions(task: TaskDefinitionElement): Seq[ExpressionElement] =
    task.commandSection.parts.flatMap(_.parts).collect { case p: PlaceholderCommandPartElement => p.expressionElement }

  private def outputExpressions(task: TaskDefinitionElement): Seq[ExpressionElement] =
    task.outputsSection.toSeq.flatMap(_.outputs).map(_.expression)

  private def otherExpressions(task: TaskDefinitionElement): Seq[ExpressionElement] =
    task.inputsSection.toSeq.flatMap(_.inputDeclarations).flatMap(_.expression) ++
      task.declarations.map(_.expression) ++
      task.runtimeSection.toSeq.flatMap(_.runtimeAttributes).map(_.value)
}
