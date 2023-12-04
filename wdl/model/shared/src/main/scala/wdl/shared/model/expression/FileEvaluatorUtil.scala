package wdl.shared.model.expression

import wom.types.WomType
import wom.values.{WomArray, WomFile, WomMap, WomObject, WomOptionalValue, WomPair, WomValue}

import scala.util.Success

object FileEvaluatorUtil {
  def findFilesToDelocalize(value: WomValue, coerceTo: WomType, coerce: Boolean = true): Seq[WomFile] = {
    val coercedValue = if (coerce) coerceTo.coerceRawValue(value) else Success(value)
    coercedValue match {
      case Success(f: WomFile) => Seq(f)
      case Success(a: WomArray) =>
        a.value.flatMap(findFilesToDelocalize(_, coerceTo, coerce = false))
      case Success(m: WomMap) =>
        (m.value flatMap { case (k, v) => Seq(k, v) } flatMap (findFilesToDelocalize(_,
                                                                                     coerceTo,
                                                                                     coerce = false
        ))).toSeq
      case Success(WomOptionalValue(_, Some(v))) => findFilesToDelocalize(v, coerceTo, coerce = false)
      case Success(WomPair(l, r)) =>
        findFilesToDelocalize(l, coerceTo, coerce = false) ++ findFilesToDelocalize(r, coerceTo, coerce = false)
      case Success(o: WomObject) =>
        o.values.values.flatMap(inner => findFilesToDelocalize(inner, inner.womType, coerce = false)).toSeq
      case _ => Seq.empty[WomFile]
    }
  }
}
