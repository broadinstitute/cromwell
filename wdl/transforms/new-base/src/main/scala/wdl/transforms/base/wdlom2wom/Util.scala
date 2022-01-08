package wdl.transforms.base.wdlom2wom

import wdl.model.draft3.elements.{MetaSectionElement, ParameterMetaSectionElement}

trait Util {

  def processMetaSections(meta: Option[MetaSectionElement], parameterMeta: Option[ParameterMetaSectionElement]) = {
    val metaMap = meta.map(_.meta).getOrElse(Map.empty)
    val parameterMetaMap = parameterMeta.map(_.metaAttributes).getOrElse(Map.empty)

    (metaMap, parameterMetaMap)
  }

}
