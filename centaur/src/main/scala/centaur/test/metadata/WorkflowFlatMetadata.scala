package centaur.test.metadata

import java.util.UUID

import cats.data.Validated._
import centaur.json.JsonUtils.EnhancedJsValue
import com.typesafe.config.Config
import common.collections.EnhancedCollections._
import common.validation.ErrorOr._
import common.validation.Validation._
import configs.Result
import configs.syntax._
import cromwell.api.model.WorkflowMetadata
import mouse.all._
import spray.json._

import scala.util.{Failure, Success, Try}

/**
  * Workflow metadata that has been flattened for Centaur test purposes. The keys are similar to the simpleton-syntax
  * stored in the Cromwell database, and values are primitive types, not nested JSON objects or arrays.
  */
case class WorkflowFlatMetadata(value: Map[String, JsValue]) extends AnyVal {

  def diff(actual: WorkflowFlatMetadata, workflowID: UUID, cacheHitUUID: Option[UUID] = None): Iterable[String] = {
    // If the test fails in initialization there wouldn't be workflow root metadata, and if that's the expectation
    // then that's ok.
    val workflowRoot = actual.value.get("workflowRoot").collectFirst { case JsString(r) => r } getOrElse "No Workflow Root"
    val missingErrors = value.keySet.diff(actual.value.keySet) map { k => s"Missing key: $k" }
    val mismatchErrors = value.keySet.intersect(actual.value.keySet) flatMap { k => diffValues(k, value(k), actual.value(k),
      workflowID, workflowRoot, cacheHitUUID)}

    mismatchErrors ++ missingErrors
  }

  private def diffValues(key: String, expected: JsValue, actual: JsValue, workflowID: UUID, workflowRoot: String, cacheHitUUID: Option[UUID]): Option[String] = {
    /*
      FIXME/TODO:

      At the moment all expected values are coerced into JsString (as it is hard to discern intended type from
      HOCON). However values coming from the metadata endpoint are JsValue. At the moment we're only using
      JsString, JsNumber and JsBoolean for comparison so this is a hacky way of handling that situation. It's
      entirely likely that it won't survive long term.
     */

    lazy val substitutedValue = expected.toString.replace("<<UUID>>", workflowID.toString).replace("<<WORKFLOW_ROOT>>", workflowRoot)
    lazy val cacheSubstitutions = cacheHitUUID match {
      case Some(uuid) => substitutedValue.replace("<<CACHE_HIT_UUID>>", uuid.toString)
      case None => substitutedValue
    }

    def stripQuotes(str: String) = str.stripPrefix("\"").stripSuffix("\"")

    val matchError = actual match {
      case o: JsString if stripQuotes(cacheSubstitutions).startsWith("~~") =>
        val stripped = stripQuotes(cacheSubstitutions).stripPrefix("~~")
        (!stripQuotes(o.toString).contains(stripped)).option(s"Actual value ${o.toString()} does not contain $stripped")
      case o: JsString => (cacheSubstitutions != o.toString).option(s"expected: $cacheSubstitutions but got: $actual")
      case o: JsNumber => (expected != JsString(o.value.toString)).option(s"expected: $cacheSubstitutions but got: $actual")
      case o: JsBoolean => (expected != JsString(o.value.toString)).option(s"expected: $cacheSubstitutions but got: $actual")
      case o: JsArray if stripQuotes(cacheSubstitutions).startsWith("~>") =>
        val stripped = stripQuotes(cacheSubstitutions).stripPrefix("~>")
        val replaced = stripped.replaceAll("\\\\\"", "\"")
        (replaced != o.toString).option(s"expected: $cacheSubstitutions but got: $actual")
      case o: JsArray => (expected != JsString(o.toString)).option(s"expected: $cacheSubstitutions but got: $actual")
      case JsNull => (expected != JsNull).option(s"expected: $cacheSubstitutions but got: $actual")
      case _ => Option(s"expected: $cacheSubstitutions but got: $actual")
    }

    matchError.map(s"Metadata mismatch for $key - " + _)
  }
}

object WorkflowFlatMetadata {
  def fromConfig(config: Config): ErrorOr[WorkflowFlatMetadata] = {
    config.extract[Map[String, Option[String]]] match {
      case Result.Success(m) => Valid(WorkflowFlatMetadata(m safeMapValues { _.map(JsString.apply).getOrElse(JsNull) }))
      case Result.Failure(_) => invalidNel(s"Metadata block can not be converted to a Map: $config")
    }
  }

  def fromMetadataJson(json: WorkflowMetadata): ErrorOr[WorkflowFlatMetadata] = {
    import DefaultJsonProtocol._
    Try(json.value.parseJson.asJsObject.flatten().convertTo[Map[String, JsValue]]) match {
      case Success(m) => Valid(WorkflowFlatMetadata(m))
      case Failure(e) => invalidNel(s"Unable to create Metadata from JSON: ${e.getMessage}")
    }
  }

  implicit class EnhancedWorkflowMetadata(val json: WorkflowMetadata) {
    def asFlat: WorkflowFlatMetadata = {
      WorkflowFlatMetadata.fromMetadataJson(json).unsafe
    }
  }
}
