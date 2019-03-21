package wdl.transforms.base.wdlom2wom

import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, _}
import wdl.model.draft3.elements.MetaSectionElement
import wom.callable.MetaValueElement

trait Util {

  def processMetaSection(maybeMetaSection: Option[MetaSectionElement]) : ErrorOr[Map[String, String]] = maybeMetaSection match {
    case None => Map.empty[String, String].validNel
    case Some(MetaSectionElement(meta)) =>
      meta traverse {
        case (key, MetaValueElement.MetaValueElementString(value)) =>
          (key -> value).validNel
        case (key, other) =>
          (key -> s"Compound types not yet supported, see #4746. String approximation: ${other.toString}").validNel
      }
  }

}
