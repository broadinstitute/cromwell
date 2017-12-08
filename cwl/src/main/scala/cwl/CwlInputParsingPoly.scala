package cwl

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import shapeless.Poly1
import wom.executable.Executable.DelayedCoercionFunction
import wom.types.{WomArrayType, WomFileType, WomType}
import wom.values._

private [cwl] object CwlInputCoercion extends Poly1 {
  implicit def cwlFileToWomValue: Case.Aux[MyriadInputValuePrimitives, DelayedCoercionFunction] = at[MyriadInputValuePrimitives] {
    _.fold(CwlInputPrimitiveCoercion)
  }
  
  implicit def inputArrayValueToWomValue: Case.Aux[Array[MyriadInputValuePrimitives], DelayedCoercionFunction] =
    at[Array[MyriadInputValuePrimitives]] { arrayValue =>
      womType: WomType => {
        import cats.instances.list._
        import cats.syntax.traverse._

        womType match {
          case womArrayType: WomArrayType =>
            arrayValue.toList
              .traverse[ErrorOr, WomValue](_.fold(CwlInputPrimitiveCoercion).apply(womArrayType.memberType))
              .map { WomArray(womArrayType, _) }

          case other => s"Cannot convert an array input value into a non array type: $other".invalidNel
        }
      }
    }
}

private [cwl] object CwlInputPrimitiveCoercion extends Poly1 {
  implicit def cwlFileToWomValue: Case.Aux[File, DelayedCoercionFunction] = at[File] { cwlFile =>
    womType: WomType => {
      womType match {
        case WomFileType => cwlFile.asWomValue
        case otherType => s"Input value is a File but the targeted input is a $otherType".invalidNel
      }
    }
  }

  implicit def stringToWomValue: Case.Aux[String, DelayedCoercionFunction] = at[String] { stringValue =>
    womType: WomType => {
      womType.coerceRawValue(stringValue).toErrorOr
    }
  }

  implicit def booleanToWomValue: Case.Aux[Boolean, DelayedCoercionFunction] = at[Boolean] { booleanValue =>
    womType: WomType => {
      womType.coerceRawValue(booleanValue).toErrorOr
    }
  }

  implicit def intToWomValue: Case.Aux[Int, DelayedCoercionFunction] = at[Int] { intValue =>
    womType: WomType => {
      womType.coerceRawValue(intValue).toErrorOr
    }
  }

  implicit def floatToWomValue: Case.Aux[Float, DelayedCoercionFunction] = at[Float] { floatValue =>
    womType: WomType => {
      womType.coerceRawValue(floatValue).toErrorOr
    }
  }

  implicit def doubleToWomValue: Case.Aux[Double, DelayedCoercionFunction] = at[Double] { doubleValue =>
    womType: WomType => {
      womType.coerceRawValue(doubleValue).toErrorOr
    }
  }

  implicit def longToWomValue: Case.Aux[Long, DelayedCoercionFunction] = at[Long] { longValue =>
    womType: WomType => {
      womType.coerceRawValue(longValue).toErrorOr
    }
  }
}
