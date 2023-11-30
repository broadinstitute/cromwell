package wdl.shared.transforms.evaluation.values

import wom.types.{WomArrayType, WomType}
import wom.values.{WomArray, WomValue}

import scala.util.{Failure, Try}

object EngineFunctions {
  def transpose(a: WomValue): Try[WomArray] = {
    case class ExpandedTwoDimensionalArray(innerType: WomType, value: Seq[Seq[WomValue]])
    def validateAndExpand(value: WomValue): Try[ExpandedTwoDimensionalArray] = value match {
      case WomArray(WomArrayType(WomArrayType(innerType)), array: Seq[WomValue]) =>
        expandWdlArray(array) map { ExpandedTwoDimensionalArray(innerType, _) }
      case WomArray(WomArrayType(nonArrayType), _) =>
        Failure(
          new IllegalArgumentException(
            s"Array must be two-dimensional to be transposed but given array of $nonArrayType"
          )
        )
      case otherValue =>
        Failure(
          new IllegalArgumentException(
            s"Function 'transpose' must be given a two-dimensional array but instead got ${otherValue.typeName}"
          )
        )
    }

    def expandWdlArray(outerArray: Seq[WomValue]): Try[Seq[Seq[WomValue]]] = Try {
      outerArray map {
        case array: WomArray => array.value
        case otherValue =>
          throw new IllegalArgumentException(
            s"Function 'transpose' must be given a two-dimensional array but instead got WdlArray[${otherValue.typeName}]"
          )
      }
    }

    def transpose(expandedTwoDimensionalArray: ExpandedTwoDimensionalArray): Try[WomArray] = Try {
      val innerType = expandedTwoDimensionalArray.innerType
      val array = expandedTwoDimensionalArray.value
      WomArray(WomArrayType(WomArrayType(innerType)), array.transpose map { WomArray(WomArrayType(innerType), _) })
    }

    validateAndExpand(a) flatMap transpose
  }
}
