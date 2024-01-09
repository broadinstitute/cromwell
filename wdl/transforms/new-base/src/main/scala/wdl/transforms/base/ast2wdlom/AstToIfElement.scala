package wdl.transforms.base.ast2wdlom

import common.validation.ErrorOr._
import cats.syntax.apply._
import cats.syntax.either._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements._

object AstToIfElement {
  def astToIfElement(implicit
    astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement],
    astNodeToWorkflowGraphElement: CheckedAtoB[GenericAstNode, WorkflowGraphElement]
  ): CheckedAtoB[GenericAst, IfElement] = CheckedAtoB.fromErrorOr("parse if block") { ast =>
    val conditionCollectionExpressionValidation: ErrorOr[ExpressionElement] =
      ast.getAttributeAs[ExpressionElement]("expression").toValidated.contextualizeErrors("parse if (...) condition")
    val bodyValidation: ErrorOr[Vector[WorkflowGraphElement]] =
      ast.getAttributeAsVector[WorkflowGraphElement]("body").toValidated.contextualizeErrors("parse 'if' { ... } body")

    (conditionCollectionExpressionValidation, bodyValidation) mapN { (condition, body) =>
      IfElement(condition, body)
    }
  }
}
