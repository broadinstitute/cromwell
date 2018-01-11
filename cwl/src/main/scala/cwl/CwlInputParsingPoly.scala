package cwl

import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Validation._
import shapeless.Poly1
import wom.executable.Executable.DelayedCoercionFunction
import wom.types._
import wom.values._

private [cwl] object CwlInputCoercion extends Poly1 {
  implicit def cwlFileToWomValue: Case.Aux[MyriadInputValuePrimitives, DelayedCoercionFunction] = at {
    _.fold(CwlInputPrimitiveCoercion)
  }

  implicit def inputArrayValueToWomValue: Case.Aux[Array[MyriadInputValuePrimitives], DelayedCoercionFunction] =
    at { arrayValue =>
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
  implicit def cwlFileOrDirectoryToWomValue: Case.Aux[FileOrDirectory, DelayedCoercionFunction] = {
    at {
      _.fold(this)
    }
  }

  implicit def cwlFileToWomValue: Case.Aux[File, DelayedCoercionFunction] = at { cwlFile =>
    womType: WomType => {
      cwlFile.asWomValue flatMap {
        womType.coerceRawValue(_).toErrorOr
      }
    }
  }

  implicit def cwlDirectoryToWomValue: Case.Aux[Directory, DelayedCoercionFunction] = at { cwlDirectory =>
    womType: WomType => {
      cwlDirectory.asWomValue flatMap {
        womType.coerceRawValue(_).toErrorOr
      }
    }
  }

  implicit def stringToWomValue: Case.Aux[String, DelayedCoercionFunction] = at { stringValue =>
    womType: WomType => {
      womType.coerceRawValue(stringValue).toErrorOr
    }
  }

  implicit def booleanToWomValue: Case.Aux[Boolean, DelayedCoercionFunction] = at { booleanValue =>
    womType: WomType => {
      womType.coerceRawValue(booleanValue).toErrorOr
    }
  }

  implicit def intToWomValue: Case.Aux[Int, DelayedCoercionFunction] = at { intValue =>
    womType: WomType => {
      womType.coerceRawValue(intValue).toErrorOr
    }
  }

  implicit def floatToWomValue: Case.Aux[Float, DelayedCoercionFunction] = at { floatValue =>
    womType: WomType => {
      womType.coerceRawValue(floatValue).toErrorOr
    }
  }

  implicit def doubleToWomValue: Case.Aux[Double, DelayedCoercionFunction] = at { doubleValue =>
    womType: WomType => {
      womType.coerceRawValue(doubleValue).toErrorOr
    }
  }

  implicit def longToWomValue: Case.Aux[Long, DelayedCoercionFunction] = at { longValue =>
    womType: WomType => {
      womType.coerceRawValue(longValue).toErrorOr
    }
  }
}
