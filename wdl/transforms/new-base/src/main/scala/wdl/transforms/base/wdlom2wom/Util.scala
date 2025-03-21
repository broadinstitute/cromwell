package wdl.transforms.base.wdlom2wom

import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.{InputsSectionElement, MetaSectionElement, OutputsSectionElement, ParameterMetaSectionElement}

trait Util {

  def processMetaSections(meta: Option[MetaSectionElement], parameterMeta: Option[ParameterMetaSectionElement]) = {
    val metaMap = meta.map(_.meta).getOrElse(Map.empty)
    val parameterMetaMap = parameterMeta.map(_.metaAttributes).getOrElse(Map.empty)

    (metaMap, parameterMetaMap)
  }

  def validateParameterMetaEntries(parameterMetaSectionElement: Option[ParameterMetaSectionElement],
                                           inputs: Option[InputsSectionElement],
                                           outputs: Option[OutputsSectionElement]
                                          ): ErrorOr[Unit] = {
    val validKeys: List[String] =
      inputs.toList.flatMap(_.inputDeclarations.map(_.name)) ++ outputs.toList.flatMap(_.outputs.map(_.name))
    val errors = parameterMetaSectionElement.toList.flatMap { pmse =>
      val keys = pmse.metaAttributes.keySet.toList
      val duplicationErrors = keys.groupBy(identity).collect {
        case (name, list) if list.size > 1 => s"Found ${list.size} parameter meta entries for '$name' (expected 0 or 1)"
      }
      val notValidKeyErrors = keys.collect {
        case name if !validKeys.contains(name) =>
          s"Invalid parameter_meta entry for '$name': not an input or output parameter"
      }
      duplicationErrors.toList ++ notValidKeyErrors
    }
    NonEmptyList.fromList(errors) match {
      case Some(nel) => Invalid(nel)
      case None => Valid(())
    }
  }

}
