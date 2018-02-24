package wdl.draft3.transforms.wdlom2wom.graph

import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, _}
import wdl.model.draft3.graph.expression.WomTypeMaker.ops._
import wdl.model.draft3.graph.expression.WomExpressionMaker.ops._
import wdl.draft3.transforms.linking.typemakers._
import wdl.draft3.transforms.linking.expression._
import wdl.model.draft3.elements._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wom.expression.WomExpression
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.ExposedExpressionNode
import wom.types.WomType

object WorkflowGraphElementToGraphNode {
  def convert(a: GraphNodeMakerInputs): ErrorOr[GraphNode] = a.node match {
    case InputDeclarationElement(typeElement, name, None) =>
      typeElement.determineWomType(Map.empty) map { womType =>
        RequiredGraphInputNode(WomIdentifier(name), womType, name)
      }
    case DeclarationElement(typeElement, name, Some(expr)) =>
      val womExprValidation: ErrorOr[WomExpression] = expr.makeWomExpression(a.linkableValues)
      val womTypeValidation: ErrorOr[WomType] = typeElement.determineWomType(Map.empty)

      (womExprValidation, womTypeValidation) flatMapN { (womExpr, womType) =>
        a.node match {
          case _: InputDeclarationElement =>
            OptionalGraphInputNodeWithDefault.apply(WomIdentifier(name), womType, womExpr, name).validNel : ErrorOr[GraphNode]
          case _: IntermediateValueDeclarationElement =>
            ExposedExpressionNode.fromInputMapping(WomIdentifier(name), womExpr, womType, a.linkablePorts)
          case _: OutputDeclarationElement =>
            ExpressionBasedGraphOutputNode.fromInputMapping(WomIdentifier(s"${a.workflowName}.$name"), womExpr, womType, a.linkablePorts)
        }

      }
  }
}

final case class GraphNodeMakerInputs(node: WorkflowGraphElement,
                                      linkableValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                                      linkablePorts: Map[String, OutputPort],
                                      workflowName: String)
