package centaur.test.metadata

import java.util.UUID

import cats.data.Validated._
import centaur.test.metadata.JsValueEnhancer._
import com.typesafe.config.Config
import common.collections.EnhancedCollections._
import common.validation.ErrorOr._
import common.validation.Validation._
import configs.Result
import configs.syntax._
import cromwell.api.model.{WorkflowId, WorkflowLabels, WorkflowMetadata, WorkflowOutputs}
import mouse.all._
import spray.json._

import scala.language.postfixOps
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
    import WorkflowFlatMetadata._

    lazy val substitutedValue = expected.toString.replaceExpectationVariables(WorkflowId(workflowID), workflowRoot)
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

  def fromWorkflowMetadata(workflowMetadata: WorkflowMetadata): ErrorOr[WorkflowFlatMetadata] = {
    val jsValue: ErrorOr[JsValue] = Try(workflowMetadata.value.parseJson) match {
      case Success(jv) => Valid(jv)
      case Failure(e) => invalidNel(s"Unable to create JsValue from String: ${e.getMessage}")
    }

    for {
      jv <- jsValue
      map <- jv.asMap
    } yield WorkflowFlatMetadata(map)
  }

  implicit class EnhancedWorkflowMetadata(val workflowMetadata: WorkflowMetadata) {
    def asFlat: WorkflowFlatMetadata = {
      WorkflowFlatMetadata.fromWorkflowMetadata(workflowMetadata).unsafe
    }
  }

  implicit class EnhancedWorkflowFlatMetadata(val workflowFlatMetadata: WorkflowFlatMetadata) {
    def stringifyValues: Map[String, JsValue] = {
      import mouse.all._
      workflowFlatMetadata.value.map {
        case (k, v: JsString) => k -> v
        case (k, JsNull) => k -> JsNull
        case (k, v) => k -> (v.toString |> JsString.apply)
      }
    }
  }

  implicit class EnhancedExpectation(val expectation: String) extends AnyVal {
    def replaceExpectationVariables(workflowId: WorkflowId, workflowRoot: String): String = {
      expectation.replaceAll("<<UUID>>", workflowId.toString).replaceAll("<<WORKFLOW_ROOT>>", workflowRoot)
    }
  }
}

object WorkflowFlatOutputs {
  implicit class EnhancedWorkflowOutputs(val workflowOutputs: WorkflowOutputs) extends AnyVal {
    def asFlat: WorkflowFlatMetadata = {
      workflowOutputs.outputs.asMap map WorkflowFlatMetadata.apply unsafe
    }
  }
}

object WorkflowFlatLabels {
  implicit class EnhancedWorkflowLabels(val workflowLabels: WorkflowLabels) extends AnyVal {
    def asFlat: WorkflowFlatMetadata = {
      workflowLabels.labels.asMap map WorkflowFlatMetadata.apply unsafe
    }
  }
}

object JsValueEnhancer {
  implicit class EnhancedJsValueAsMap(val jsValue: JsValue) extends AnyVal {
    import DefaultJsonProtocol._
    import centaur.json.JsonUtils._

    def asMap: ErrorOr[Map[String, JsValue]] = {
      Try(jsValue.asJsObject.flatten().convertTo[Map[String, JsValue]]) match {
        case Success(m) => Valid(m)
        case Failure(e) => invalidNel(s"Unable to convert JsValue to JsObject: ${e.getMessage}")
      }
    }
  }
}
