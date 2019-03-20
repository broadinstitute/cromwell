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
          s"non string meta values are not handled currently <${key} -> ${other}>, see https://github.com/broadinstitute/cromwell/issues/4746".invalidNel
      }
  }

}
