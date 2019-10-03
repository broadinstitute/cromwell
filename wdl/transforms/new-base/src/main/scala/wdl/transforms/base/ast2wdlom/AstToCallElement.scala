package wdl.transforms.base.ast2wdlom

import cats.syntax.validated._
import cats.syntax.either._
import cats.instances.either._
import cats.syntax.apply._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import wdl.model.draft3.elements.{CallBodyElement, CallElement}
import wdl.model.draft3.elements.ExpressionElement.KvPair
import wom.SourceFileLocation
object AstToCallElement {

  def astToCallElement(implicit astNodeToKvPair: CheckedAtoB[GenericAstNode, KvPair]): CheckedAtoB[GenericAst, CallElement] = CheckedAtoB.fromErrorOr { ast =>
    def convertBodyElement(a: GenericAst): Checked[CallBodyElement] = {
      a.getAttributeAsVector[KvPair]("inputs") map CallBodyElement
    }

    val callableNameValidation: ErrorOr[String] = astNodeToString(ast.getAttribute("task")).toValidated

    val aliasValidation: ErrorOr[Option[String]] = Option(ast.getAttribute("alias")) match {
      case Some(a) => astNodeToString(a).map(Some.apply).toValidated
      case None => None.validNel
    }

    val afterValidation: ErrorOr[Vector[String]] = ast.getAttributeAsVector[String]("after", optional = true).toValidated

    implicit val astNodeToCallBodyElement: CheckedAtoB[GenericAstNode, CallBodyElement] = astNodeToAst andThen CheckedAtoB.fromCheck(convertBodyElement _)

    val callBodyValidation: ErrorOr[Option[CallBodyElement]] = ast.getAttributeAsOptional[CallBodyElement]("body").toValidated

    val sourceLocation : Option[SourceFileLocation] = ast.getSourceLine.map(SourceFileLocation(_))

    // This 'mapN' is split into two so that if we have a call name we can include it in the error message
    (callableNameValidation, aliasValidation) flatMapN { (name, alias) =>
      val result = (afterValidation, callBodyValidation) mapN { (after, body) =>
        CallElement(name, alias, after, body, sourceLocation)
      }
      result.contextualizeErrors(s"call $name" + alias.fold("")(a => s" as $a"))
    }
  }

}
