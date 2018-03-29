package wdl.draft3.transforms.wdlom2wom.graph

import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, _}
import wdl.model.draft3.graph.expression.WomTypeMaker.ops._
import wdl.model.draft3.graph.expression.WomExpressionMaker.ops._
import wdl.draft3.transforms.linking.typemakers._
import wdl.draft3.transforms.linking.expression._
import wdl.model.draft3.elements._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedCallOutputOrIdentifierAndMemberAccessHook, UnlinkedConsumedValueHook, UnlinkedIdentifierHook}
import wdl.shared.transforms.wdlom2wom.WomGraphMakerTools
import wom.callable.{CallableTaskDefinition, Callable, CommandTaskDefinition, TaskDefinition}
import wom.expression.WomExpression
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.ExposedExpressionNode
import wom.types.WomType

object WorkflowGraphElementToGraphNode {
  def convert(a: GraphNodeMakerInputs): ErrorOr[Set[GraphNode]] = a.node match {
    case ie: InputDeclarationElement =>
      val inputNodeMakerInputs = GraphInputNodeMakerInputs(ie, a.linkableValues, a.linkablePorts, a.availableTypeAliases, a.workflowName)
      InputDeclarationElementToGraphNode.convert(inputNodeMakerInputs)

    case DeclarationElement(typeElement, name, Some(expr)) =>
      val womExprValidation: ErrorOr[WomExpression] = expr.makeWomExpression(a.availableTypeAliases, a.linkableValues)
      val womTypeValidation: ErrorOr[WomType] = typeElement.determineWomType(a.availableTypeAliases)

      (womExprValidation, womTypeValidation) flatMapN { (womExpr, womType) =>
        a.node match {
          case _: InputDeclarationElement =>
            Set[GraphNode](OptionalGraphInputNodeWithDefault.apply(WomIdentifier(name), womType, womExpr, name)).validNel
          case _: IntermediateValueDeclarationElement =>
            ExposedExpressionNode.fromInputMapping(WomIdentifier(name), womExpr, womType, a.linkablePorts) map { Set(_) }
          case _: OutputDeclarationElement =>
            ExpressionBasedGraphOutputNode.fromInputMapping(WomIdentifier(name, s"${a.workflowName}.$name"), womExpr, womType, a.linkablePorts) map {Set(_)}
        }
      }

    case se: ScatterElement =>
      val scatterMakerInputs = ScatterNodeMakerInputs(se, a.linkableValues, a.linkablePorts, a.availableTypeAliases, a.workflowName, a.insideAScatter, a.callables)
      ScatterElementToGraphNode.convert(scatterMakerInputs)

    case ie: IfElement =>
      val ifMakerInputs = ConditionalNodeMakerInputs(ie, a.linkableValues, a.linkablePorts, a.availableTypeAliases, a.workflowName, a.insideAScatter, a.callables)
      IfElementToGraphNode.convert(ifMakerInputs)

    case ce: CallElement =>
      val callableNodeMakerInputs = CallableNodeMakerInputs(ce, a.linkableValues, a.linkablePorts, a.availableTypeAliases, a.workflowName, a.insideAScatter, a.callables)
      CallElementToGraphNode.convert(callableNodeMakerInputs)
  }
}

final case class GraphNodeMakerInputs(node: WorkflowGraphElement,
                                      linkableValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                                      linkablePorts: Map[String, OutputPort],
                                      availableTypeAliases: Map[String, WomType],
                                      workflowName: String,
                                      insideAScatter: Boolean,
                                      callables: Set[Callable])
