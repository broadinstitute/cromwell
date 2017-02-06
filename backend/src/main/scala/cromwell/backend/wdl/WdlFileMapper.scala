package cromwell.backend.wdl

import lenthall.util.TryUtil
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlOptionalValue, WdlPair, WdlValue}

import scala.util.{Success, Try}

object WdlFileMapper {
  def mapWdlFiles(mapper: (WdlFile => WdlFile))(wdlValue: WdlValue): Try[WdlValue] = {
    wdlValue match {
      case file: WdlFile => Try(mapper(file))
      case array: WdlArray =>
        val mappedArray = array.value map mapWdlFiles(mapper)
        TryUtil.sequence(mappedArray) map {
          WdlArray(array.wdlType, _)
        }
      case map: WdlMap =>
        val mappedMap = map.value map {
          case (key, value) => mapWdlFiles(mapper)(key) -> mapWdlFiles(mapper)(value)
        }

        TryUtil.sequenceKeyValues(mappedMap) map {
          WdlMap(map.wdlType, _)
        }
      case pair: WdlPair =>
        val mappedPair: (Try[WdlValue], Try[WdlValue]) = (mapWdlFiles(mapper)(pair.left), mapWdlFiles(mapper)(pair.right))
        TryUtil.sequenceTuple(mappedPair) map {
          (WdlPair.apply _).tupled
        }
      case optionalValue: WdlOptionalValue =>
        val mappedOptional: Option[Try[WdlValue]] = optionalValue.value.map(mapWdlFiles(mapper))
        TryUtil.sequenceOption(mappedOptional) map {
          WdlOptionalValue(optionalValue.innerType, _)
        }
      case other => Success(other)
    }
  }
}
