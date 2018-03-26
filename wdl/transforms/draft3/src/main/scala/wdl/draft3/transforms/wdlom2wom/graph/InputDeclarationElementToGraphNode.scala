package wdl.draft3.transforms.wdlom2wom.graph

import cats.syntax.apply._
import common.validation.ErrorOr.{ErrorOr, _}
import wdl.model.draft3.graph.expression.WomTypeMaker.ops._
import wdl.model.draft3.graph.expression.WomExpressionMaker.ops._
import wdl.draft3.transforms.linking.typemakers._
import wdl.draft3.transforms.linking.expression._
import wdl.draft3.transforms.linking.expression.consumed._
import wdl.model.draft3.elements.ExpressionElement.{ArrayLiteral, IdentifierLookup, SelectFirst}
import wdl.model.draft3.elements._
import wdl.model.draft3.graph._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wom.expression.WomExpression
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.ExposedExpressionNode
import wom.types.{WomOptionalType, WomType}

object InputDeclarationElementToGraphNode {
  def convert(a: GraphInputNodeMakerInputs): ErrorOr[Set[GraphNode]] = a.node match {
    case InputDeclarationElement(typeElement, name, None) =>
      val nameInInputSet = s"${a.workflowName}.$name"
      typeElement.determineWomType(a.availableTypeAliases) map {
        case opt: WomOptionalType => Set(OptionalGraphInputNode(WomIdentifier(name), opt.flatOptionalType, nameInInputSet))
        case womType => Set(RequiredGraphInputNode(WomIdentifier(name), womType, nameInInputSet))
      }
    case InputDeclarationElement(typeElement, name, Some(expr)) if expr.expressionConsumedValueHooks.isEmpty =>
      val womExprValidation: ErrorOr[WomExpression] = expr.makeWomExpression(a.availableTypeAliases, a.linkableValues)
      val womTypeValidation: ErrorOr[WomType] = typeElement.determineWomType(a.availableTypeAliases)
      (womExprValidation, womTypeValidation) mapN { (womExpr, womType) =>
        Set[GraphNode](OptionalGraphInputNodeWithDefault.apply(WomIdentifier(name), womType, womExpr, name))
      }
    // This input must have upstream dependencies but WOM won't allow that, so make two WOM nodes:
    // - An optional input graph node
    // - An expression node representing select_first([optional_input, original_expression])
    case InputDeclarationElement(typeElement, name, Some(expr)) =>
      val womTypeValidation: ErrorOr[WomType] = typeElement.determineWomType(a.availableTypeAliases)

      val nameInInputSet = s"${a.workflowName}.$name"

      womTypeValidation flatMap { womType =>
        val optionalInputIdentifier = s"__$name"
        val optionalInputType = WomOptionalType(womType).flatOptionalType
        val optionalInputNode = OptionalGraphInputNode(WomIdentifier(optionalInputIdentifier), optionalInputType, nameInInputSet)
        val optionalInputNodeGeneratedValue = GeneratedIdentifierValueHandle(optionalInputIdentifier, optionalInputType)
        val collectorExpressionElement: ExpressionElement = SelectFirst(ArrayLiteral(Seq(IdentifierLookup(optionalInputIdentifier), expr)))
        val collectorExpressionExtraConsumedValueEntry = UnlinkedIdentifierHook(optionalInputIdentifier) -> optionalInputNodeGeneratedValue
        val collectorLinkablePorts = a.linkablePorts + (optionalInputIdentifier -> optionalInputNode.singleOutputPort)

        val collectorWomExpressionValidation: ErrorOr[WomExpression] = collectorExpressionElement.makeWomExpression(a.availableTypeAliases, a.linkableValues + collectorExpressionExtraConsumedValueEntry)
        val collectorNode: ErrorOr[GraphNode] = collectorWomExpressionValidation flatMap { womExpr =>
          ExposedExpressionNode.fromInputMapping(WomIdentifier(name, s"${a.workflowName}.$name"), womExpr, womType, collectorLinkablePorts)
        }

        collectorNode map {
          Set(optionalInputNode, _)
        }
      }
  }
}

final case class GraphInputNodeMakerInputs(node: InputDeclarationElement,
                                            linkableValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                                            linkablePorts: Map[String, OutputPort],
                                            availableTypeAliases: Map[String, WomType],
                                            workflowName: String)
