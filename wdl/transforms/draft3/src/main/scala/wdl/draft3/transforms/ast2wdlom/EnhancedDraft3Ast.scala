package wdl.draft3.transforms.ast2wdlom

import cats.data.NonEmptyList
import cats.instances.vector._
import cats.instances.either._
import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import wdl.draft3.parser.WdlParser.{Ast, AstList, AstNode}

import scala.collection.JavaConverters._

object EnhancedDraft3Ast {

  implicit class EnhancedAstNode(val astNode: AstNode) extends AnyVal {
    def astListAsVector: Checked[Vector[AstNode]] = astNode match {
      case list: AstList => list.asScala.toVector.validNelCheck
      case _ => s"Invalid target for astListAsVector: ${astNode.getClass.getSimpleName}".invalidNelCheck
    }
  }

  implicit class EnhancedAst(val ast: Ast) extends AnyVal {

    def getAttributeAsAstNodeVector(attr: String): Checked[Vector[AstNode]] = for {
      attributeNode <- Option(ast.getAttribute(attr)).toChecked(s"No attribute $attr found on Ast of type ${ast.getName}")
      asVector <- attributeNode.astListAsVector
    } yield asVector

    /**
      * Will get an attribute on this Ast as an AstList and then convert that into a vector of Ast
      * @param attr The attribute to read from this Ast
      */
    def getAttributeAsVector[A](attr: String)(implicit toA: CheckedAtoB[AstNode, A]): Checked[Vector[A]] = {
      for {
        asVector <- getAttributeAsAstNodeVector(attr)
        result <- asVector.traverse(toA.run)
      } yield result
    }

    def getAttributeAsVectors[A, B](attr: String, astName1: String, astName2: String)
                                   (implicit astNodeToA: CheckedAtoB[AstNode, A], astNodeToB: CheckedAtoB[AstNode, B]
                                   ): (Checked[(Vector[A], Vector[B])]) = {

      ast.getAttributeAsVector[Ast](attr) flatMap { asts =>

        val validated = ( collectByAstName(asts, astName1).traverse[ErrorOr, A] { type1Ast => astNodeToA.run(type1Ast).toValidated },
          collectByAstName(asts, astName2).traverse[ErrorOr, B] { type2Asts => astNodeToB.run(type2Asts).toValidated },
          badAstsValidation(asts, Set(astName1, astName2))
        ) mapN {
          (type1Values, type2Values, _) => (type1Values, type2Values)
        }
        validated.toEither
      }
    }

    private def collectByAstName(asts: Vector[Ast], astName: String) = asts collect { case yes if yes.getName == astName => yes }
    private def badAstsValidation(asts: Vector[Ast], validAstNames: Set[String]): ErrorOr[Unit] = {
      def badAstDescription(ast: Ast) = s"Unexpected AST type: ${ast.getName} (expected one of: ${validAstNames.mkString(", ")})"
      asts collect {
        case no if !validAstNames.contains(no.getName) => no
      } match {
        case empty if empty.isEmpty => ().validNel
        case nonEmpty => (NonEmptyList.fromListUnsafe(nonEmpty.toList) map badAstDescription).invalid
      }
    }
  }
}
