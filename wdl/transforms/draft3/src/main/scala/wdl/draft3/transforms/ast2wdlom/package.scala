package wdl.draft3.transforms

import cats.instances.vector._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import wdl.draft3.parser.WdlParser.{Ast, AstList, AstNode}

import scala.collection.JavaConverters._

package object ast2wdlom {

  implicit val draft3FileElementFromAst = Draft3FileElementFromAst
  implicit val draft3ImportElementFromAst = Draft3ImportElementFromAst
  implicit val draft3TaskDefinitionElementFromAst = Draft3TaskDefinitionElementFromAst
  implicit val draft3WorkflowDefinitionElementFromAst = Draft3WorkflowDefinitionElementFromAst

  implicit class EnhancedAstNode(val astNode: AstNode) extends AnyVal {
    def astListAsVector: Checked[Vector[AstNode]] = astNode match {
      case list: AstList => list.asScala.toVector.validNelCheck
      case _ => s"Invalid target for astListAsVector: ${astNode.getClass.getSimpleName}".invalidNelCheck
    }
  }

  implicit class EnhancedAst(val ast: Ast) extends AnyVal {

    /**
      * Will get an attribute on this Ast as an AstList and then convert that into a vector of Ast
      * @param attr The attribute to read from this Ast
      */
    def getAttributeAsAstVector(attr: String): Checked[Vector[Ast]] = {

      val astNodesCheck = for {
        attributeNode <- Option(ast.getAttribute(attr)).toChecked(s"No attribute $attr found on Ast of type ${ast.getName}")
        asVector <- attributeNode.astListAsVector
      } yield asVector

      def astNodeToAst(astNode: AstNode): ErrorOr[Ast] = astNode match {
        case a: Ast => a.validNel
        case other => s"Cannot cast type '${other.getClass.getSimpleName}' into Ast (while getting attribute $attr on ${ast.getName})".invalidNel
      }

      astNodesCheck flatMap { _.traverse(astNodeToAst).toEither }
    }
  }
}
