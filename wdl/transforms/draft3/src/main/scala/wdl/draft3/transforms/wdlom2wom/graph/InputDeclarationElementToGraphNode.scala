package wdl.draft3.transforms.wdlom2wom.graph

import cats.syntax.apply._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.graph.expression.WomTypeMaker.ops._
import wdl.model.draft3.graph.expression.WomExpressionMaker.ops._
import wdl.draft3.transforms.linking.typemakers._
import wdl.draft3.transforms.linking.expression._
import wdl.model.draft3.elements._
import wdl.model.draft3.graph._
import wom.expression.WomExpression
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.types.{WomOptionalType, WomType}

object InputDeclarationElementToGraphNode {
  def convert(a: GraphInputNodeMakerInputs): ErrorOr[Set[GraphNode]] = a.node match {
    case InputDeclarationElement(typeElement, name, None) =>
      val nameInInputSet = s"${a.workflowName}.$name"
      typeElement.determineWomType(a.availableTypeAliases) map {
        case opt: WomOptionalType => Set(OptionalGraphInputNode(WomIdentifier(name), opt.flatOptionalType, nameInInputSet))
        case womType => Set(RequiredGraphInputNode(WomIdentifier(name), womType, nameInInputSet))
      }
    case InputDeclarationElement(typeElement, name, Some(expr)) =>
      // NOTE: This won't work if the input expression has upstream dependencies. That's why we had to eliminate them
      // in the WorkflowDefinitionElementToWomWorkflowDefinition convert step.
      val womExprValidation: ErrorOr[WomExpression] = expr.makeWomExpression(a.availableTypeAliases, a.linkableValues)
      val womTypeValidation: ErrorOr[WomType] = typeElement.determineWomType(a.availableTypeAliases)
      (womExprValidation, womTypeValidation) mapN { (womExpr, womType) =>
        val nameInInputSet = s"${a.workflowName}.$name"
        Set[GraphNode](OptionalGraphInputNodeWithDefault.apply(WomIdentifier(name), womType, womExpr, nameInInputSet))
      }
  }
}

final case class GraphInputNodeMakerInputs(node: InputDeclarationElement,
                                            linkableValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                                            linkablePorts: Map[String, OutputPort],
                                            availableTypeAliases: Map[String, WomType],
                                            workflowName: String)
