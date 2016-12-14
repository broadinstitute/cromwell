package cromwell.backend.validation

import cats.data.ValidatedNel
import cats.syntax.validated._
import cromwell.core.{EvaluatedRuntimeAttributes, OptionNotFoundException, WorkflowOptions}
import lenthall.util.TryUtil
import wdl4s.types.WdlType
import wdl4s.values.WdlValue

import scala.util.{Failure, Try}

object RuntimeAttributesDefault {

  def workflowOptionsDefault(options: WorkflowOptions, mapping: Map[String, Traversable[WdlType]]):
  Try[Map[String, WdlValue]] = {
    options.defaultRuntimeOptions flatMap { attrs =>
      TryUtil.sequenceMap(attrs collect {
        case (k, v) if mapping.contains(k) =>
          val maybeTriedValue = mapping(k) map {  _.coerceRawValue(v) } find { _.isSuccess } getOrElse {
            Failure(new RuntimeException(s"Could not parse JsonValue $v to valid WdlValue for runtime attribute $k"))
          }
          k -> maybeTriedValue
      }, "Failed to coerce default runtime options")
    } recover {
      case _: OptionNotFoundException => Map.empty[String, WdlValue]
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
