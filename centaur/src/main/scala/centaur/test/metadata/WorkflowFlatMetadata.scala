package centaur.test.metadata

import java.util.UUID
import cats.data.Validated._
import centaur.test.metadata.JsValueEnhancer._
import com.typesafe.config.{Config, ConfigValue, ConfigValueType}
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.api.model.{WorkflowId, WorkflowLabels, WorkflowMetadata, WorkflowOutputs}
import mouse.all._
import spray.json._

import java.util
import scala.jdk.CollectionConverters._
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
    val workflowRoot =
      actual.value.get("workflowRoot").collectFirst { case JsString(r) => r } getOrElse "No Workflow Root"
    val missingErrors = value.keySet.diff(actual.value.keySet) map { k => s"Missing key: $k" }
    val mismatchErrors = value.keySet.intersect(actual.value.keySet) flatMap { k =>
      diffValues(k, value(k), actual.value(k), workflowID, workflowRoot, cacheHitUUID)
    }

    mismatchErrors ++ missingErrors
  }

  private def diffValues(key: String,
                         expected: JsValue,
                         actual: JsValue,
                         workflowID: UUID,
                         workflowRoot: String,
                         cacheHitUUID: Option[UUID]
  ): Option[String] = {
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
      case o: JsString if stripQuotes(cacheSubstitutions).endsWith("~~") =>
        val stripped = stripQuotes(cacheSubstitutions).stripSuffix("~~")
        (!stripQuotes(o.toString).contains(stripped))
          .option(s"Actual value ${o.toString()} does not start with $stripped")
      case o: JsString =>
        (cacheSubstitutions != o.toString).option(s"expected: $cacheSubstitutions but got: $actual")
      case o: JsNumber =>
        expected match {
          case JsNumber(value) =>
            (value.compare(o.value) != 0).option(s"expected: $cacheSubstitutions but got: $actual")
          case _ => Option(s"metadata $actual is a number, but expected to be the same type as $expected")
        }
      case o: JsBoolean =>
        expected match {
          case JsBoolean(value) => (value != o.value).option(s"expected: $cacheSubstitutions but got: $actual")
          case _ => Option(s"metadata $actual is a boolean, but expected to be the same type as $expected")
        }
      case o: JsArray if stripQuotes(cacheSubstitutions).startsWith("~>") =>
        val stripped = stripQuotes(cacheSubstitutions).stripPrefix("~>")
        val replaced = stripped.replaceAll("\\\\\"", "\"")
        replaced.parseJson match {
          case JsArray(elements) =>
            (elements != o.elements).option(s"expected: $replaced.parseJson but got: $actual")
          case _ => Option(s"metadata $actual is an array, but expected to be the same type as ${replaced.parseJson}")
        }
      case o: JsArray =>
        expected match {
          case JsArray(elements) =>
            (elements != o.elements).option(s"expected: $cacheSubstitutions but got: $actual")
          case _ => Option(s"metadata $actual is an array, but expected to be the same type as $expected")
        }
      case JsNull => (expected != JsNull).option(s"expected: $cacheSubstitutions but got: $actual")
      case _ => Option(s"expected: $cacheSubstitutions but got: $actual")
    }

    matchError.map(s"Metadata mismatch for $key - " + _)
  }
}

object WorkflowFlatMetadata {

  def fromConfig(config: Config): ErrorOr[WorkflowFlatMetadata] = {
    val metadataValues = config.entrySet().toArray(Array[util.Map.Entry[String, ConfigValue]]()).map { value =>
      val v: ConfigValue = value.getValue
      val key: String = value.getKey.stripPrefix("\"").stripSuffix("\"")
      v.valueType() match {
        case ConfigValueType.BOOLEAN => key -> JsBoolean.apply(v.unwrapped().asInstanceOf[Boolean])
        case ConfigValueType.NUMBER =>
          val num = v.unwrapped().asInstanceOf[Number]
          if (num.intValue() == num.doubleValue()) {
            key -> JsNumber.apply(num.intValue())
          } else {
            key -> JsNumber.apply(num.doubleValue())
          }
        case ConfigValueType.LIST =>
          key ->
            parseArray(
              v.unwrapped()
                .asInstanceOf[util.List[Object]]
                .asScala
                .toVector
            )
        case ConfigValueType.NULL => key -> JsNull
        case _ =>
          /*
            FIXME/TODO:

            As part of WX-1629, this was changed to allow Centaur tests to perform comparisons that take data type
            into account instead of casting all values to strings before comparing. However, it seems that expected
            values from *.test files are read in as either Strings, Numbers, or Lists -- as such, in order to compare
            booleans, a cast is required. This isn't ideal because it makes it impossible to use the Strings "true" or
            "false" -- they will always be cast to booleans. The changes necessary to fix this are more wide reaching
            than the scope of WX-1629, so this must be fixed in the future.
           */
          if (v.unwrapped().asInstanceOf[String] == "true" || v.unwrapped().asInstanceOf[String] == "false") {
            key -> JsBoolean.apply(v.unwrapped().asInstanceOf[String].toBoolean)
          } else {
            key -> JsString.apply(v.unwrapped().asInstanceOf[String])
          }
      }
    }

    Valid(WorkflowFlatMetadata(metadataValues.toMap))
  }

  private def parseArray(element: Any): JsValue =
    element match {
      case array: Vector[Any] =>
        JsArray(
          array.map {
            case value: java.util.ArrayList[_] => parseArray(value.asInstanceOf[util.List[Any]].asScala.toVector)
            case value: Any => parseArray(value)
          }
        )
      case _ =>
        element match {
          case value: java.lang.Boolean => JsBoolean.apply(value)
          case value: java.lang.Integer => JsNumber.apply(value)
          case value: java.lang.Float => JsNumber.apply(value)
          case value: java.lang.Double => JsNumber.apply(value)
          case _ => JsString.apply(element.asInstanceOf[String])
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
    def asFlat: WorkflowFlatMetadata =
      WorkflowFlatMetadata.fromWorkflowMetadata(workflowMetadata).unsafe
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
    def replaceExpectationVariables(workflowId: WorkflowId, workflowRoot: String): String =
      expectation.replaceAll("<<UUID>>", workflowId.toString).replaceAll("<<WORKFLOW_ROOT>>", workflowRoot)
  }
}

object WorkflowFlatOutputs {
  implicit class EnhancedWorkflowOutputs(val workflowOutputs: WorkflowOutputs) extends AnyVal {
    def asFlat: WorkflowFlatMetadata =
      workflowOutputs.outputs.asMap map WorkflowFlatMetadata.apply unsafe
  }
}

object WorkflowFlatLabels {
  implicit class EnhancedWorkflowLabels(val workflowLabels: WorkflowLabels) extends AnyVal {
    def asFlat: WorkflowFlatMetadata =
      workflowLabels.labels.asMap map WorkflowFlatMetadata.apply unsafe
  }
}

object JsValueEnhancer {
  implicit class EnhancedJsValueAsMap(val jsValue: JsValue) extends AnyVal {
    import DefaultJsonProtocol._
    import centaur.json.JsonUtils._

    def asMap: ErrorOr[Map[String, JsValue]] =
      Try(jsValue.asJsObject.flatten().convertTo[Map[String, JsValue]]) match {
        case Success(m) => Valid(m)
        case Failure(e) => invalidNel(s"Unable to convert JsValue to JsObject: ${e.getMessage}")
      }
  }
}
