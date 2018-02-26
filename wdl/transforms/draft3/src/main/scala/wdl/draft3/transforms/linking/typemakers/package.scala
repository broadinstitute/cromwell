package wdl.draft3.transforms.linking

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.{PrimitiveTypeElement, TypeElement}
import wdl.model.draft3.graph.expression.WomTypeMaker
import wdl.model.draft3.graph.expression.WomTypeMaker.ops._
import wom.types.WomType

package object typemakers {
  implicit object primitiveTypeElementConverter extends WomTypeMaker[PrimitiveTypeElement] {
    override def determineWomType(a: PrimitiveTypeElement, availableAliases: Map[String, WomType]): ErrorOr[WomType] = {
      a.primitiveType.validNel
    }
  }

  implicit object typeElementToWomType extends WomTypeMaker[TypeElement] {
    override def determineWomType(a: TypeElement, availableAliases: Map[String, WomType]): ErrorOr[WomType] = a match {
      case p: PrimitiveTypeElement => p.determineWomType(availableAliases)
      case other => s"No rule to convert type element '$other' to WOM type".invalidNel
    }
  }
}
