package centaur.test.metadata

import java.util.UUID

import centaur.test.ErrorOr
import cats.data.Validated._
import com.typesafe.config.Config
import configs.syntax._
import spray.json._
import centaur.json.JsonUtils.EnhancedJsValue
import configs.Result

import scala.util.{Failure, Success, Try}

case class WorkflowMetadata(value: Map[String, JsValue]) extends AnyVal {

  def diff(other: WorkflowMetadata, workflowID: UUID): Iterable[String] = {
    val missingErrors = value.keySet.diff(other.value.keySet) map { k => s"Missing key: $k" }
    val mismatchErrors = value.keySet.intersect(other.value.keySet) flatMap { k => diffValues(k, value(k), other.value(k), workflowID) }

    mismatchErrors ++ missingErrors
  }

  private def diffValues(key: String, expected: JsValue, other: JsValue, workflowID: UUID): Option[String] = {
    /*
      FIXME/TODO:

      At the moment all expected values are coerced into JsString (as it is hard to discern intended type from
      HOCON). However values coming from the metadata endpoint are JsValue. At the moment we're only using
      JsString, JsNumber and JsBoolean for comparison so this is a hacky way of handling that situation. It's
      entirely likely that it won't survive long term.
     */

    lazy val sanitisedExpectedUUID = expected.toString.replace("<<UUID>>", workflowID.toString)

    val isMatch = other match {
      case o: JsString => sanitisedExpectedUUID == o.toString
      case o: JsNumber => expected == JsString(o.value.toString)
      case o: JsBoolean => expected == JsString(o.value.toString)
      case _ => false
    }

    if (isMatch) None
    else Option(s"Metadata mismatch for $key - expected: $sanitisedExpectedUUID but got: $other")
  }

}

object WorkflowMetadata {
  def fromConfig(config: Config): ErrorOr[WorkflowMetadata] = {
    config.extract[Map[String, String]] match {
      case Result.Success(m) => Valid(WorkflowMetadata(m mapValues { JsString(_) }))
      case Result.Failure(_) => invalidNel("Metadata block can not be converted to a Map")
    }
  }

  def fromMetadataJson(json: String): ErrorOr[WorkflowMetadata] = {
    import DefaultJsonProtocol._

    Try(json.parseJson.asJsObject.flatten().convertTo[Map[String, JsValue]]) match {
      case Success(m) => Valid(WorkflowMetadata(m))
      case Failure(e) => invalidNel(s"Unable to create Metadata from JSON: ${e.getMessage}")
    }
  }

}
