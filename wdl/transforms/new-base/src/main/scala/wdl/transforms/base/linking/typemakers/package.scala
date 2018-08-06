package wdl.transforms.base.linking

import cats.syntax.apply._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import common.validation.Validation._
import wdl.model.draft3.elements._
import wdl.model.draft3.graph.expression.WomTypeMaker
import wdl.model.draft3.graph.expression.WomTypeMaker.ops._
import wom.types._

package object typemakers {
  implicit val primitiveTypeElementConverter: WomTypeMaker[PrimitiveTypeElement] = new WomTypeMaker[PrimitiveTypeElement] {
    override def determineWomType(a: PrimitiveTypeElement, availableAliases: Map[String, WomType]): ErrorOr[WomType] = {
      a.primitiveType.validNel
    }
  }

  implicit val arrayTypeElementConverter: WomTypeMaker[ArrayTypeElement] = new WomTypeMaker[ArrayTypeElement] {
    override def determineWomType(a: ArrayTypeElement, availableAliases: Map[String, WomType]): ErrorOr[WomType] = {
      a.inner.determineWomType(availableAliases) map { inner => WomArrayType(inner) }
    }
  }

  implicit val mapTypeElementConverter: WomTypeMaker[MapTypeElement] = new WomTypeMaker[MapTypeElement] {
    override def determineWomType(a: MapTypeElement, availableAliases: Map[String, WomType]): ErrorOr[WomType] = {
      (a.keyType.determineWomType(availableAliases),
        a.valueType.determineWomType(availableAliases)) mapN { (keyType, valueType) => WomMapType(keyType, valueType) }
    }
  }

  implicit val optionalTypeElementConverter: WomTypeMaker[OptionalTypeElement] = new WomTypeMaker[OptionalTypeElement] {
    override def determineWomType(a: OptionalTypeElement, availableAliases: Map[String, WomType]): ErrorOr[WomType] = {
      a.maybeType.determineWomType(availableAliases) map { inner => WomOptionalType(inner) }
    }
  }

  implicit val nonEmptyTypeElementConverter: WomTypeMaker[NonEmptyTypeElement] = new WomTypeMaker[NonEmptyTypeElement] {
    override def determineWomType(a: NonEmptyTypeElement, availableAliases: Map[String, WomType]): ErrorOr[WomType] = {
      a.arrayType.determineWomType(availableAliases) flatMap {
        case WomArrayType(memberType) => WomNonEmptyArrayType(memberType).validNel
        case other: WomType => s"Cannot declare a non-empty $other (+ is only applicable to Array[_] types)".invalidNel
      }
    }
  }

  implicit val pairTypeElementConverter: WomTypeMaker[PairTypeElement] = new WomTypeMaker[PairTypeElement] {
    override def determineWomType(a: PairTypeElement, availableAliases: Map[String, WomType]): ErrorOr[WomType] = {
      (a.leftType.determineWomType(availableAliases),
        a.rightType.determineWomType(availableAliases)) mapN { (keyType, valueType) => WomPairType(keyType, valueType) }
    }
  }

  implicit val structTypeElementConverter: WomTypeMaker[TypeAliasElement] = new WomTypeMaker[TypeAliasElement] {
    override def determineWomType(a: TypeAliasElement, availableAliases: Map[String, WomType]): ErrorOr[WomType] = {
      availableAliases.get(a.alias).toErrorOr(s"No struct definition for '${a.alias}' found in available structs: [${availableAliases.values.mkString(", ")}]")
    }
  }

  implicit val typeElementToWomType: WomTypeMaker[TypeElement] = new WomTypeMaker[TypeElement] {
    override def determineWomType(a: TypeElement, availableAliases: Map[String, WomType]): ErrorOr[WomType] = a match {
      case p: PrimitiveTypeElement => p.determineWomType(availableAliases)
      case a: ArrayTypeElement => a.determineWomType(availableAliases)
      case m: MapTypeElement => m.determineWomType(availableAliases)
      case o: OptionalTypeElement => o.determineWomType(availableAliases)
      case n: NonEmptyTypeElement => n.determineWomType(availableAliases)
      case p: PairTypeElement => p.determineWomType(availableAliases)
      case s: TypeAliasElement => s.determineWomType(availableAliases)
      case ObjectTypeElement => WomObjectType.validNel

      case other => s"No rule to convert type element '$other' to WOM type".invalidNel
    }
  }
}
