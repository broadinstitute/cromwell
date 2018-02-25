package cromwell.backend.validation

import cats.data.ValidatedNel
import cats.syntax.validated._
import cromwell.core.{OptionNotFoundException, WorkflowOptions}
import common.util.TryUtil
import wom.core.EvaluatedRuntimeAttributes
import wom.types.WomType
import wom.values.WomValue

import scala.util.{Failure, Try}

object RuntimeAttributesDefault {

  def workflowOptionsDefault(options: WorkflowOptions, mapping: Map[String, Traversable[WomType]]):
  Try[Map[String, WomValue]] = {
    options.defaultRuntimeOptions flatMap { attrs =>
      TryUtil.sequenceMap(attrs collect {
        case (k, v) if mapping.contains(k) =>
          val maybeTriedValue = mapping(k) map {  _.coerceRawValue(v) } find { _.isSuccess } getOrElse {
            Failure(new RuntimeException(s"Could not parse JsonValue $v to valid WomValue for runtime attribute $k"))
          }
          k -> maybeTriedValue
      }, "Failed to coerce default runtime options")
    } recover {
      case _: OptionNotFoundException => Map.empty[String, WomValue]
    }
  }

  /**
    * Traverse defaultsList in order, and for each of them add the missing (and only missing) runtime attributes.
   */
  def withDefaults(attrs: EvaluatedRuntimeAttributes, defaultsList: List[EvaluatedRuntimeAttributes]): EvaluatedRuntimeAttributes = {
    defaultsList.foldLeft(attrs)((acc, default) => {
      acc ++ default.filterKeys(!acc.keySet.contains(_))
    })
  }

  def noValueFoundFor[A](attribute: String): ValidatedNel[String, A] = s"Can't find an attribute value for key $attribute".invalidNel
}
