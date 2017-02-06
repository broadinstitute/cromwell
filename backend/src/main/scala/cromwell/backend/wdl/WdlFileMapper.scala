package cromwell.backend.wdl

import lenthall.util.TryUtil
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlValue}

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
          // TODO: TryUtil.sequenceMap doesn't currently handle keys of Try[WdlValue], so for now calling .get here
          case (key, value) => mapWdlFiles(mapper)(key).get -> mapWdlFiles(mapper)(value)
        }
        TryUtil.sequenceMap(mappedMap) map {
          WdlMap(map.wdlType, _)
        }
      case other => Success(other)
    }
  }
}
