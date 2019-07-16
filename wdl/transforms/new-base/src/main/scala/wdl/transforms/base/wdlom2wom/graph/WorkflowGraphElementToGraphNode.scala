package wdl.transforms.base.wdlom2wom.graph

import cats.syntax.apply._
import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, _}
import wdl.model.draft3.graph.expression.WomTypeMaker.ops._
import wdl.model.draft3.graph.expression.WomExpressionMaker.ops._
import wdl.transforms.base.linking.typemakers._
import wdl.transforms.base.linking.expression._
import wdl.model.draft3.elements._
import wdl.model.draft3.graph.expression.{FileEvaluator, TypeEvaluator, ValueEvaluator}
import wdl.model.draft3.graph.{ExpressionValueConsumer, GeneratedValueHandle, UnlinkedConsumedValueHook}
import wom.callable.Callable
import wom.expression.WomExpression
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.ExposedExpressionNode
import wom.types.WomType
import wdl.transforms.base.wdlom2wdl.WdlWriter.ops._
import wdl.transforms.base.wdlom2wdl.WdlWriterImpl._

object WorkflowGraphElementToGraphNode {
  def convert(a: GraphNodeMakerInputs)
             (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement],
              fileEvaluator: FileEvaluator[ExpressionElement],
              typeEvaluator: TypeEvaluator[ExpressionElement],
              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[GraphNode]] = a.node match {
    case ie: InputDeclarationElement =>
      val inputNodeMakerInputs = GraphInputNodeMakerInputs(ie, a.linkableValues, a.linkablePorts, a.availableTypeAliases, a.workflowName)
      InputDeclarationElementToGraphNode.convert(inputNodeMakerInputs)

    case DeclarationElement(typeElement, name, Some(expr)) =>
      val womExprValidation: ErrorOr[WomExpression] = expr.makeWomExpression(a.availableTypeAliases, a.linkableValues)
      val womTypeValidation: ErrorOr[WomType] = typeElement.determineWomType(a.availableTypeAliases)

      val result = (womExprValidation, womTypeValidation) flatMapN { (womExpr, womType) =>
        val correctType: ErrorOr[Unit] = validateAssignmentType(womExpr, womType)

        val graphNode: ErrorOr[Set[GraphNode]] = a.node match {
          case _: InputDeclarationElement =>
            Set[GraphNode](OptionalGraphInputNodeWithDefault.apply(WomIdentifier(name), womType, womExpr, name)).validNel
          case _: IntermediateValueDeclarationElement =>
            ExposedExpressionNode.fromInputMapping(WomIdentifier(name), womExpr, womType, a.linkablePorts) map { Set(_) }
          case _: OutputDeclarationElement =>
            ExpressionBasedGraphOutputNode.fromInputMapping(WomIdentifier(name, s"${a.workflowName}.$name"), womExpr, womType, a.linkablePorts) map {Set(_)}
        }

        (correctType, graphNode) mapN { (_, gn) => gn }
      }
      result.contextualizeErrors(s"process declaration '${typeElement.toWdlV1} $name = ${expr.toWdlV1}'")

    case se: ScatterElement =>
      val scatterMakerInputs = ScatterNodeMakerInputs(se, a.upstreamCalls, a.linkableValues, a.linkablePorts, a.availableTypeAliases, a.workflowName, a.insideAScatter, a.convertNestedScatterToSubworkflow, a.callables)
      ScatterElementToGraphNode.convert(scatterMakerInputs)

    case ie: IfElement =>
      val ifMakerInputs = ConditionalNodeMakerInputs(ie, a.upstreamCalls, a.linkableValues, a.linkablePorts, a.availableTypeAliases, a.workflowName, a.insideAScatter, a.convertNestedScatterToSubworkflow, a.callables)
      IfElementToGraphNode.convert(ifMakerInputs)

    case ce: CallElement =>
      val callNodeMakerInputs = CallNodeMakerInputs(ce, a.upstreamCalls, a.linkableValues, a.linkablePorts, a.availableTypeAliases, a.workflowName, a.insideAScatter, a.callables)
      CallElementToGraphNode.convert(callNodeMakerInputs)
  }

  def validateAssignmentType(womExpr: WomExpression, womType: WomType): ErrorOr[Unit] = {
    // NB: WdlomWomExpressions don't actually use the inputTypes argument, they already know their type:
    val exprType: ErrorOr[WomType] = womExpr.evaluateType(inputTypes = Map.empty)
    val correctType: ErrorOr[Unit] = exprType.flatMap {
      case goodType if womType.isCoerceableFrom(goodType) => ().validNel
      case badType => s"Cannot coerce expression of type '${badType.stableName}' to '${womType.stableName}'".invalidNel
    }
    correctType
  }
}

final case class GraphNodeMakerInputs(node: WorkflowGraphElement,
                                      upstreamCalls: Map[String, CallNode],
                                      linkableValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                                      linkablePorts: Map[String, OutputPort],
                                      availableTypeAliases: Map[String, WomType],
                                      workflowName: String,
                                      insideAScatter: Boolean,
                                      convertNestedScatterToSubworkflow : Boolean,
                                      callables: Map[String, Callable])
