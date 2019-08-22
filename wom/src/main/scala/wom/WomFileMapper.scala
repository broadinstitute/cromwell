package wom

import java.io.FileNotFoundException

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
  def mapWomFiles(mapper: WomFile => WomFile,
                  exceptions: Set[WomFile])
                 (womValue: WomValue): Try[WomValue] = {
    womValue match {
      case file: WomFile if !exceptions.contains(file) => Try(mapper(file))
      case array: WomArray =>
        val mappedArray = array.value map mapWomFiles(mapper, exceptions)
        TryUtil.sequence(mappedArray) map {
          WomArray(array.womType, _)
        }
      case map: WomMap =>
        val mappedMap = map.value map {
          case (key, value) => mapWomFiles(mapper, exceptions)(key) -> mapWomFiles(mapper, exceptions)(value)
        }
        TryUtil.sequenceKeyValues(mappedMap) map {
          WomMap(map.womType, _)
        }
      case womObject: WomObjectLike =>
        val mappedMap = womObject.values map {
          case (key, value) => key -> mapWomFiles(mapper, exceptions)(value)
        }
        TryUtil.sequenceMap(mappedMap).map(WomObject.withTypeUnsafe(_, womObject.womObjectTypeLike))
      case pair: WomPair =>
        val mappedPair: (Try[WomValue], Try[WomValue]) = (mapWomFiles(mapper, exceptions)(pair.left), mapWomFiles(mapper, exceptions)(pair.right))
        TryUtil.sequenceTuple(mappedPair) map {
          (WomPair.apply _).tupled
        }
      case optionalValue: WomOptionalValue =>
        // Build a `WomOptionalValue` from an `Option[WomValue]`.
        def buildWomOptionalValue(optionalWomValue: Option[WomValue]) = WomOptionalValue(optionalValue.innerType, optionalWomValue)

        val mappedOptional: Option[Try[WomValue]] = optionalValue.value.map(mapWomFiles(mapper, exceptions))
        mappedOptional match {
          case Some(o) =>
            o map Option.apply recover { case _: FileNotFoundException => None } map buildWomOptionalValue
          case None => Success(buildWomOptionalValue(None))
        }
      case coproduct: WomCoproductValue => mapWomFiles(mapper, exceptions)(coproduct.womValue)
      case other => Success(other)
    }
  }
}
