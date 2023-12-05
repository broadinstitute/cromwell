package wdl.transforms.base.ast2wdlom

import cats.syntax.either._
import cats.data.NonEmptyList
import common.transforms.CheckedAtoB
import common.validation.Checked._
import wdl.model.draft3.elements._

object AstToWorkflowGraphNodeElementConverterMaker {
  private def astToWorkflowGraphNodeElement(
    astNodeToDeclarationContent: CheckedAtoB[GenericAstNode, DeclarationContent],
    astNodeToCallElement: CheckedAtoB[GenericAstNode, CallElement],
    astNodeToScatterElement: CheckedAtoB[GenericAstNode, ScatterElement],
    astNodeToIfElement: CheckedAtoB[GenericAstNode, IfElement]
  ): CheckedAtoB[GenericAst, WorkflowGraphElement] = CheckedAtoB.fromCheck { a: GenericAst =>
    a.getName match {
      case "Declaration" => astNodeToDeclarationContent(a).map(IntermediateValueDeclarationElement.fromContent)
      case "Call" => astNodeToCallElement(a)
      case "Scatter" => astNodeToScatterElement(a)
      case "If" => astNodeToIfElement(a)
      case other => s"No conversion defined for Ast with name $other to WorkflowGraphElement".invalidNelCheck
    }
  }
}

class AstToWorkflowGraphNodeElementConverterMaker() {
  var astNodeToDeclarationContent: Option[CheckedAtoB[GenericAstNode, DeclarationContent]] = None
  var astNodeToCallElement: Option[CheckedAtoB[GenericAstNode, CallElement]] = None
  var astNodeToScatterElement: Option[CheckedAtoB[GenericAstNode, ScatterElement]] = None
  var astNodeToIfElement: Option[CheckedAtoB[GenericAstNode, IfElement]] = None

  val uninitializedMessage = NonEmptyList.fromListUnsafe(List(""))

  implicit val astNodeToDeclarationContentFixed: CheckedAtoB[GenericAstNode, DeclarationContent] =
    CheckedAtoB.fromCheck { a =>
      Either.fromOption(astNodeToDeclarationContent, uninitializedMessage) flatMap { c => c.run(a) }
    }

  implicit val astNodeToCallElementFixed: CheckedAtoB[GenericAstNode, CallElement] = CheckedAtoB.fromCheck { a =>
    Either.fromOption(astNodeToCallElement, uninitializedMessage) flatMap { c => c.run(a) }
  }

  implicit val astNodeToScatterElementFixed: CheckedAtoB[GenericAstNode, ScatterElement] = CheckedAtoB.fromCheck { a =>
    Either.fromOption(astNodeToScatterElement, uninitializedMessage) flatMap { c => c.run(a) }
  }

  implicit val astNodeToIfElementFixed: CheckedAtoB[GenericAstNode, IfElement] = CheckedAtoB.fromCheck { a =>
    Either.fromOption(astNodeToIfElement, uninitializedMessage) flatMap { c => c.run(a) }
  }

  val converter = AstToWorkflowGraphNodeElementConverterMaker.astToWorkflowGraphNodeElement(
    astNodeToDeclarationContentFixed,
    astNodeToCallElementFixed,
    astNodeToScatterElementFixed,
    astNodeToIfElementFixed
  )
}
