package cwl

import cats.syntax.validated._
import lenthall.validation.ErrorOr.ErrorOr
import lenthall.validation.Validation._
import shapeless.Poly1
import wom.executable.Executable.DelayedCoercionFunction
import wom.types.{WdlArrayType, WdlFileType, WdlType}
import wom.values._

private [cwl] object CwlInputCoercion extends Poly1 {
  implicit def cwlFileToWdlValue: Case.Aux[MyriadInputValuePrimitives, DelayedCoercionFunction] = at[MyriadInputValuePrimitives] {
    _.fold(CwlInputPrimitiveCoercion)
  }
  
  implicit def inputArrayValueToWdlValue: Case.Aux[Array[MyriadInputValuePrimitives], DelayedCoercionFunction] =
    at[Array[MyriadInputValuePrimitives]] { arrayValue =>
      womType: WdlType => {
        import cats.instances.list._
        import cats.syntax.traverse._

        womType match {
          case wdlArrayType: WdlArrayType =>
            arrayValue.toList
              .traverse[ErrorOr, WdlValue](_.fold(CwlInputPrimitiveCoercion).apply(wdlArrayType.memberType))
              .map { WdlArray(wdlArrayType, _) }

          case other => s"Cannot convert an array input value into a non array type: $other".invalidNel
        }
      }
    }
}

private [cwl] object CwlInputPrimitiveCoercion extends Poly1 {
  implicit def cwlFileToWdlValue: Case.Aux[File, DelayedCoercionFunction] = at[File] { cwlFile =>
    womType: WdlType => {
      womType match {
        case WdlFileType => cwlFile.asWdlValue
        case otherType => s"Input value is a File but the targeted input is a $otherType".invalidNel
      }
    }
  }

  implicit def stringToWdlValue: Case.Aux[String, DelayedCoercionFunction] = at[String] { stringValue =>
    womType: WdlType => {
      womType.coerceRawValue(stringValue).toErrorOr
    }
  }

  implicit def booleanToWdlValue: Case.Aux[Boolean, DelayedCoercionFunction] = at[Boolean] { booleanValue =>
    womType: WdlType => {
      womType.coerceRawValue(booleanValue).toErrorOr
    }
  }

  implicit def intToWdlValue: Case.Aux[Int, DelayedCoercionFunction] = at[Int] { intValue =>
    womType: WdlType => {
      womType.coerceRawValue(intValue).toErrorOr
    }
  }

  implicit def floatToWdlValue: Case.Aux[Float, DelayedCoercionFunction] = at[Float] { floatValue =>
    womType: WdlType => {
      womType.coerceRawValue(floatValue).toErrorOr
    }
  }

  implicit def doubleToWdlValue: Case.Aux[Double, DelayedCoercionFunction] = at[Double] { doubleValue =>
    womType: WdlType => {
      womType.coerceRawValue(doubleValue).toErrorOr
    }
  }

  implicit def longToWdlValue: Case.Aux[Long, DelayedCoercionFunction] = at[Long] { longValue =>
    womType: WdlType => {
      womType.coerceRawValue(longValue).toErrorOr
    }
  }
}
