package wom

import common.util.TryUtil
import wom.values._

import scala.util.{Success, Try}

object WomFileMapper {
  /**
    * Loops over a WomValue applying the supplied mapper function whenever a WomFile is encountered.
    *
    * Where WomValue.collectAsSeq applies a partial function that collects into a Seq, this method applies a function to
    * all encountered WomFile within a WomMap, WomPair, WomArray, etc. Similarly WomFile.mapFile only traverses
    * the passed in WomFile plus any additional files referenced within the WomFile.
    *
    * @see [[wom.values.WomValue.collectAsSeq]]
    * @see [[wom.values.WomFile.mapFile]]
    */
  def mapWomFiles(mapper: (WomFile => WomFile))(womValue: WomValue): Try[WomValue] = {
    womValue match {
      case file: WomFile => Try(mapper(file))
      case array: WomArray =>
        val mappedArray = array.value map mapWomFiles(mapper)
        TryUtil.sequence(mappedArray) map {
          WomArray(array.womType, _)
        }
      case map: WomMap =>
        val mappedMap = map.value map {
          case (key, value) => mapWomFiles(mapper)(key) -> mapWomFiles(mapper)(value)
        }
        TryUtil.sequenceKeyValues(mappedMap) map {
          WomMap(map.womType, _)
        }
      case pair: WomPair =>
        val mappedPair: (Try[WomValue], Try[WomValue]) = (mapWomFiles(mapper)(pair.left), mapWomFiles(mapper)(pair.right))
        TryUtil.sequenceTuple(mappedPair) map {
          (WomPair.apply _).tupled
        }
      case optionalValue: WomOptionalValue =>
        val mappedOptional: Option[Try[WomValue]] = optionalValue.value.map(mapWomFiles(mapper))
        TryUtil.sequenceOption(mappedOptional) map {
          WomOptionalValue(optionalValue.innerType, _)
        }
      case other => Success(other)
    }
  }
}
