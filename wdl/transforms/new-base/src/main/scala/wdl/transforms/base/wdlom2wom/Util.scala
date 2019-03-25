package wdl.transforms.base.wdlom2wom

import wdl.model.draft3.elements.{MetaSectionElement, ParameterMetaSectionElement}
import wom.callable.MetaValueElement

trait Util {

  def processMetaSections(meta: Option[MetaSectionElement], parameterMeta: Option[ParameterMetaSectionElement]) = {

    def stringifyMetaValues(meta: Map[String, MetaValueElement]): Map[String, String] = {
      meta map {
        case (key, MetaValueElement.MetaValueElementString(value)) =>
          key -> value
        case (key, other) =>
          key -> s"Compound types not yet supported, see #4746. String approximation: ${other.toString}"
      }
    }

    val metaMap = meta.map(_.meta).getOrElse(Map.empty)
    val parameterMetaMap = parameterMeta.map(_.metaAttributes).getOrElse(Map.empty)

    (stringifyMetaValues(metaMap), stringifyMetaValues(parameterMetaMap))
  }

}
