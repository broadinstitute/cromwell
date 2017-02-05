package cromwell.backend.wdl

import lenthall.exception.AggregatedException
import lenthall.util.TryUtil
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlOptionalValue, WdlPair, WdlValue}

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
        EnhancedTryUtil.sequenceMapKeyValues(mappedMap) map {
          WdlMap(map.wdlType, _)
        }
      case pair: WdlPair =>
        val mappedPair: (Try[WdlValue], Try[WdlValue]) = (mapWdlFiles(mapper)(pair.left), mapWdlFiles(mapper)(pair.right))
        EnhancedTryUtil.sequenceMapTuple(mappedPair) map {
          (WdlPair.apply _).tupled
        }
      case optionalValue: WdlOptionalValue =>
        val mappedOptional: Option[Try[WdlValue]] = optionalValue.value.map(mapWdlFiles(mapper))
        EnhancedTryUtil.sequenceOptionalValue(mappedOptional) map {
          WdlOptionalValue(optionalValue.innerType, _)
        }
      case other => Success(other)
    }
  }
}

// TODO: Move into Lenthall
object EnhancedTryUtil {
  private def sequenceIterable[T](tries: Iterable[Try[_]], unbox: () => T, prefixErrorMessage: String) = {
    tries collect { case f: Failure[_] => f } match {
      case failures if failures.nonEmpty =>
        val exceptions = failures.toSeq.map(_.exception)
        Failure(AggregatedException(prefixErrorMessage, exceptions.toList))
      case _ => Success(unbox())
    }
  }

  def sequenceOptionalValue[T](tried: Option[Try[T]], prefixErrorMessage: String = ""): Try[Option[T]] = {
    def unbox: Option[T] = tried.map(_.get)

    val triesSeq: Seq[Try[_]] = tried.toSeq
    sequenceIterable(triesSeq, unbox _, prefixErrorMessage)
  }

  def sequenceMapTuple[T, U](tries: (Try[T], Try[U]), prefixErrorMessage: String = ""): Try[(T, U)] = {
    def unbox: (T, U) = tries match {
      case (try1, try2) => (try1.get, try2.get)
    }

    val triesSeq: Seq[Try[_]] = tries match {
      case (try1, try2) => Seq(try1, try2)
    }
    sequenceIterable(triesSeq, unbox _, prefixErrorMessage)
  }

  // TODO: TryUtil.sequenceMap doesn't currently handle keys of Try[WdlValue], so for now implementing our own here.
  def sequenceMapKeyValues[T, U](tries: Map[Try[T], Try[U]], prefixErrorMessage: String = ""): Try[Map[T, U]] = {
    def unbox: Map[T, U] = tries map { case (tryKey, tryValue) => tryKey.get -> tryValue.get }

    val triesSeq: Seq[Try[_]] = tries.toSeq.flatMap(Function.tupled(Seq(_, _)))
    sequenceIterable(triesSeq, unbox _, prefixErrorMessage)
  }
}
