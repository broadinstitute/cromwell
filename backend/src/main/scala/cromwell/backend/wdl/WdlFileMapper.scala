package cromwell.backend.wdl

import common.util.TryUtil
import wom.values._

import scala.util.{Success, Try}

object WdlFileMapper {
  def mapWdlFiles(mapper: (WomFile => WomFile))(womValue: WomValue): Try[WomValue] = {
    womValue match {
      case file: WomFile => Try(mapper(file))
      case array: WomArray =>
        val mappedArray = array.value map mapWdlFiles(mapper)
        TryUtil.sequence(mappedArray) map {
          WomArray(array.womType, _)
        }
      case map: WomMap =>
        val mappedMap = map.value map {
          case (key, value) => mapWdlFiles(mapper)(key) -> mapWdlFiles(mapper)(value)
        }
        TryUtil.sequenceKeyValues(mappedMap) map {
          WomMap(map.womType, _)
        }
      case pair: WomPair =>
        val mappedPair: (Try[WomValue], Try[WomValue]) = (mapWdlFiles(mapper)(pair.left), mapWdlFiles(mapper)(pair.right))
        TryUtil.sequenceTuple(mappedPair) map {
          (WomPair.apply _).tupled
        }
      case optionalValue: WomOptionalValue =>
        val mappedOptional: Option[Try[WomValue]] = optionalValue.value.map(mapWdlFiles(mapper))
        TryUtil.sequenceOption(mappedOptional) map {
          WomOptionalValue(optionalValue.innerType, _)
        }
      case other => Success(other)
    }
  }
}
