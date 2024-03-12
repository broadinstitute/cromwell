package wdl.transforms.base.wdlom2wom.graph

import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.graph.expression.WomTypeMaker.ops._
import wdl.model.draft3.graph.expression.WomExpressionMaker.ops._
import wdl.transforms.base.linking.typemakers._
import wdl.transforms.base.linking.expression._
import wdl.transforms.base.wdlom2wom.graph.WorkflowGraphElementToGraphNode.validateAssignmentType
import wdl.model.draft3.elements._
import wdl.model.draft3.graph._
import wdl.model.draft3.graph.expression.{FileEvaluator, TypeEvaluator, ValueEvaluator}
import wom.expression.WomExpression
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.types.{WomOptionalType, WomType}
import wdl.transforms.base.wdlom2wdl.WdlWriter.ops._
import wdl.transforms.base.wdlom2wdl.WdlWriterImpl.{expressionElementWriter, typeElementWriter}

object InputDeclarationElementToGraphNode {
  def convert(a: GraphInputNodeMakerInputs)(implicit
    expressionValueConsumer: ExpressionValueConsumer[ExpressionElement],
    fileEvaluator: FileEvaluator[ExpressionElement],
    typeEvaluator: TypeEvaluator[ExpressionElement],
    valueEvaluator: ValueEvaluator[ExpressionElement]
  ): ErrorOr[Set[GraphNode]] = a.node match {
    case InputDeclarationElement(typeElement, name, None) =>
      val nameInInputSet = s"${a.workflowName}.$name"
      typeElement.determineWomType(a.availableTypeAliases) map {
        case opt: WomOptionalType =>
          Set(OptionalGraphInputNode(WomIdentifier(name), opt.flatOptionalType, nameInInputSet))
        case womType => Set(RequiredGraphInputNode(WomIdentifier(name), womType, nameInInputSet))
      }
    case InputDeclarationElement(typeElement, name, Some(expr)) =>
      // NOTE: This won't work if the input expression has upstream dependencies. That's why we had to eliminate them
      // in the WorkflowDefinitionElementToWomWorkflowDefinition convert step.
      val womExprValidation: ErrorOr[WomExpression] = expr.makeWomExpression(a.availableTypeAliases, a.linkableValues)
      val womTypeValidation: ErrorOr[WomType] = typeElement.determineWomType(a.availableTypeAliases)
      val result = (womExprValidation, womTypeValidation) flatMapN { (womExpr, womType) =>
        validateAssignmentType(womExpr, womType) map { _ =>
          val nameInInputSet = s"${a.workflowName}.$name"
          Set[GraphNode](OptionalGraphInputNodeWithDefault.apply(WomIdentifier(name), womType, womExpr, nameInInputSet))
        }
      }
      result.contextualizeErrors(s"process input declaration '${typeElement.toWdlV1} $name = ${expr.toWdlV1}'")
  }
}

final case class GraphInputNodeMakerInputs(node: InputDeclarationElement,
                                           linkableValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                                           linkablePorts: Map[String, OutputPort],
                                           availableTypeAliases: Map[String, WomType],
                                           workflowName: String
)
