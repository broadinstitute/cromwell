package wdl.transforms.base.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.validated._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements._
import wom.SourceFileLocation

object AstToScatterElement {
  def astToScatterElement(implicit astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement],
                          astNodeToWorkflowGraphElement: CheckedAtoB[GenericAstNode, WorkflowGraphElement]
                         ): CheckedAtoB[GenericAst, ScatterElement] = CheckedAtoB.fromErrorOr("read scatter section") { ast =>

    val scatterVariableValidation: ErrorOr[GenericTerminal] = ast.getAttributeAs[GenericTerminal]("item").toValidated

    val scatterCollectionExpressionValidation: ErrorOr[ExpressionElement] = ast.getAttributeAs[ExpressionElement]("collection").toValidated
    val bodyValidation: ErrorOr[Vector[WorkflowGraphElement]] = ast.getAttributeAsVector[WorkflowGraphElement]("body").toValidated
    val srcLocValidation : ErrorOr[Option[SourceFileLocation]] = ast.getSourceLine.map(SourceFileLocation.convert).validNel

    (scatterVariableValidation, scatterCollectionExpressionValidation, bodyValidation, srcLocValidation) mapN { (variable, collection, body, srcLoc) =>
      val scatterName = s"ScatterAt${variable.getLine}_${variable.getColumn}"
      ScatterElement(scatterName, collection, variable.getSourceString, body, srcLoc)
    }
  }
}
