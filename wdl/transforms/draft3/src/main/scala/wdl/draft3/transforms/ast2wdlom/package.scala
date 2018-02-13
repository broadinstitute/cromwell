package wdl.draft3.transforms

import better.files.File
import cats.data.Kleisli
import cats.instances.either._
import common.Checked
import common.validation.Checked._
import wdl.draft3.parser.WdlParser.{Ast, AstNode, Terminal}
import wdl.draft3.transforms.ast2wdlom.CheckedAstNodeToAst.CheckedAstNodeToAst
import wdl.draft3.transforms.parsing.FileParser
import wdl.model.draft3.elements.FileElement

package object ast2wdlom {

  type CheckedAtoB[A, B] = Kleisli[Checked, A, B]
  object CheckedAtoB {
    def apply[A, B](implicit runner: CheckedAtoB[A, B]): Kleisli[Checked, A, B] = runner
    def apply[A, B](run: A => Checked[B]): Kleisli[Checked, A, B] = Kleisli(run)
  }

  type CheckedAstNodeTo[B] = CheckedAtoB[AstNode, B]
  object CheckedAstNodeTo {
    def apply[B](run: AstNode => Checked[B]): CheckedAstNodeTo[B] = CheckedAtoB[AstNode, B](run)
  }

  type CheckedAstTo[B] = CheckedAtoB[Ast, B]
  object CheckedAstTo {
    def apply[B](run: Ast => Checked[B]): CheckedAstTo[B] = CheckedAtoB[Ast, B](run)
  }

  implicit val astFromAstNode: CheckedAstNodeToAst = CheckedAstNodeToAst.instance
  implicit val draft3FileElementFromAstNode = astFromAstNode andThen CheckedAstToFileElement.instance
  implicit val draft3ImportElementFromAstNode = astFromAstNode andThen CheckedAstToImportElement.instance
  implicit val draft3TaskDefinitionElementFromAstNode = astFromAstNode andThen CheckedAstToTaskDefinitionElement.instance
  implicit val draft3WorkflowDefinitionElementFromAstNode = astFromAstNode andThen CheckedAstToWorkflowDefinitionElement.instance

  implicit val draft3FileElementFromFile: CheckedAtoB[File, FileElement] = FileParser.instance andThen CheckedAstToFileElement.instance

  implicit val checkedAstNodeToString: CheckedAstNodeTo[String] = CheckedAtoB[AstNode, String] (run = (a: AstNode) => a match {
    case t: Terminal => t.getSourceString.validNelCheck
    case other: AstNode => s"Cannot convert ${other.getClass.getSimpleName} into String".invalidNelCheck
  })

}
