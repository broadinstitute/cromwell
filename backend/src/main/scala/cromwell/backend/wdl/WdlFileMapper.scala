package cromwell.backend.wdl

import lenthall.exception.AggregatedException
import lenthall.util.TryUtil
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlValue}

import scala.util.{Failure, Success, Try}

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
        KeyValueTryUtil.sequenceMapKeyValues(mappedMap) map {
          WdlMap(map.wdlType, _)
        }
      case other => Success(other)
    }
  }
}

object KeyValueTryUtil {
  private def sequenceIterable[T](tries: Iterable[Try[_]], unbox: () => T, prefixErrorMessage: String) = {
    tries collect { case f: Failure[_] => f } match {
      case failures if failures.nonEmpty =>
        val exceptions = failures.toSeq.map(_.exception)
        Failure(AggregatedException(prefixErrorMessage, exceptions.toList))
      case _ => Success(unbox())
    }
  }

  // TODO: TryUtil.sequenceMap doesn't currently handle keys of Try[WdlValue], so for now implementing our own here.
  def sequenceMapKeyValues[T, U](tries: Map[Try[T], Try[U]], prefixErrorMessage: String = ""): Try[Map[T, U]] = {
    def unbox: Map[T, U] = tries map { case (tryKey, tryValue) => tryKey.get -> tryValue.get }

    val triesSeq: Seq[Try[_]] = tries.toSeq.flatMap(Function.tupled(Seq(_, _)))
    sequenceIterable(triesSeq, unbox _, prefixErrorMessage)
  }
}
