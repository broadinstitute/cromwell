package wdl.transforms.base.wdlom2wom

import wom.callable.MetaValueElement

trait Util {

  def processMetaSection(meta: Map[String, MetaValueElement]): Map[String, String] = {
    meta map {
      case (key, MetaValueElement.MetaValueElementString(value)) =>
        key -> value
      case (key, other) =>
        key -> s"Compound types not yet supported, see #4746. String approximation: ${other.toString}"
    }
  }

}
