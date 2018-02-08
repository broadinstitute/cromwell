package wdl.draft3.transforms.ast2wdlom

import cats.data.NonEmptyList
import cats.instances.vector._
import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import wdl.draft3.parser.WdlParser.Ast

trait FromAst[A] extends FromAtoB[Ast, A]

object FromAst {

  def apply[A](ast: Ast)(implicit value: FromAtoB[Ast, A]): ErrorOr[A] = { value.convert(ast) }

  implicit class GetAttributeAsVector(ast: Ast) {
    def getAttributeAsVector[A](attr: String)(implicit v: FromAtoB[Ast, A]): ErrorOr[Vector[A]] = {
      val astVectorValidation = ast.getAttributeAsAstVector(attr).toValidated
      astVectorValidation.flatMap { _.traverse[ErrorOr, A] { im => FromAst[A](im) } }
    }

    def getAttributeAsVectors[A, B](attr: String, astName1: String, astName2: String)
                                   (implicit astToA: FromAtoB[Ast, A], astToB: FromAtoB[Ast, B]
                                   ): (ErrorOr[(Vector[A], Vector[B])]) = {

      ast.getAttributeAsAstVector(attr).toValidated flatMap { asts =>

        ( collectByAstName(asts, astName1).traverse[ErrorOr, A] { type1Ast => FromAst[A](type1Ast) },
          collectByAstName(asts, astName2).traverse[ErrorOr, B] { type2Asts => FromAst[B](type2Asts) },
          badAstsValidation(asts, Set(astName1, astName2))
        ) mapN {
          (type1Values, type2Values, _) => (type1Values, type2Values)
        }
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

