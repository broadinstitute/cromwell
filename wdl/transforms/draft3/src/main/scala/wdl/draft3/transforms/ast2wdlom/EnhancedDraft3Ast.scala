package wdl.draft3.transforms.ast2wdlom

import cats.data.NonEmptyList
import cats.instances.map
import cats.instances.vector._
import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.flatMap
import cats.syntax.traverse._
import cats.syntax.validated._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import wdl.draft3.parser.WdlParser.{Ast, AstList, AstNode}

import scala.collection._
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
      attributeNode <- Option(ast.getAttribute(attr)).toChecked(s"No attribute '$attr' found on Ast of type ${ast.getName}. Did you mean: ${ast.getAttributes.asScala.keys.mkString(", ")}")
      asVector <- attributeNode.astListAsVector
    } yield asVector

    /**
      * Will get an attribute on this Ast as an AstNode and then convert that into a single element of
      * the required type.
      */
    def getAttributeAs[A](attr: String)(implicit toA: CheckedAtoB[AstNode, A]): Checked[A] = {
      val attribute = Option(ast.getAttribute(attr))
      attribute.map(toA.run).getOrElse(s"No attribute '$attr' found on Ast '${ast.getName}'. Did you mean: ${ast.getAttributes.asScala.keys.mkString(", ")}".invalidNelCheck)
    }

    /**
      * Will get an attribute on this Ast as an AstList and then convert that into a vector of Ast
      * @param attr The attribute to read from this Ast
      */
    def getAttributeAsVector[A](attr: String)(implicit toA: CheckedAtoB[AstNode, A]): Checked[Vector[A]] = {
      for {
        asVector <- getAttributeAsAstNodeVector(attr)
        // This toValidated/toEither dance is necessary to
        // (1) collect all errors from the traverse as an ErrorOr, then
        // (2) convert back into a Checked for the flatMap
        result <- asVector.traverse[ErrorOr, A](item => toA.run(item).toValidated).toEither
      } yield result
    }

    /**
      * Gets an attribute on this Ast as an Optional Ast, returns an empty Option if the attribute is empty.
      */
    def getAttributeAsOptional[A](attr: String)(implicit toA: CheckedAtoB[AstNode, Option[A]]): Checked[Option[A]] =  {
      val attribute: Option[AstNode] = Option(ast.getAttribute(attr))
      attribute.map(toA.run).getOrElse(Option.empty.validNelCheck)
    }
  }
}
