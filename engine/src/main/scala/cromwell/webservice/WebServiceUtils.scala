package cromwell.webservice

import akka.http.scaladsl.marshalling.ToEntityMarshaller
import akka.http.scaladsl.model._
import akka.http.scaladsl.model.headers.RawHeader
import akka.http.scaladsl.server.Directives.complete
import akka.http.scaladsl.server.Route
import akka.stream.Materializer
import akka.util.{ByteString, Timeout}
import cromwell.webservice.routes.CromwellApiService._
import io.circe._
import io.circe.generic.auto._
import io.circe.syntax._

import scala.concurrent.{ExecutionContext, Future}

trait WebServiceUtils {

  type MaterializedFormData = Map[String, ByteString]

  def materializeFormData(formData: Multipart.FormData)(implicit timeout: Timeout, materializer: Materializer, executionContext: ExecutionContext): Future[MaterializedFormData] = {
    formData.parts.mapAsync[(String, ByteString)](1) {
      bodyPart => bodyPart.toStrict(timeout.duration)(materializer).map(strict => bodyPart.name -> strict.entity.data)(executionContext)
    }.runFold(Map.empty[String, ByteString])((map, tuple) => map + tuple)(materializer)
  }

  /**
    * Completes a response of a Product, probably a case class, using an implicit marshaller, probably a json encoder.
    */
  def completeResponse[A <: Product](statusCode: StatusCode, value: A, warnings: Seq[String])
                                    (implicit mt: ToEntityMarshaller[A]): Route = {
    complete((statusCode, warningHeaders(warnings), value))
  }

  /**
    * Completes a response of string with the supplied content type.
    *
    * This is currently only used for pretty printing json, which should ideally be a query string option wired into
    * the json marshaller. See cromwell.webservice.WebServiceUtils.EnhancedThrowable below.
    *
    * Ex: https://www.elastic.co/guide/en/elasticsearch/reference/6.6/common-options.html#_pretty_results
    */
  def completeResponse(statusCode: StatusCode,
                       contentType: ContentType.NonBinary,
                       value: String,
                       warnings: Seq[String]): Route = {
    complete((statusCode, warningHeaders(warnings), HttpEntity(contentType, value)))
  }

  def warningHeaders(warnings: Seq[String]): List[HttpHeader] = {
    warnings.toList map { warning =>
      /*
      Need a quoted string.
      https://stackoverflow.com/questions/7886782

      Using a poor version of ~~#!
      https://github.com/akka/akka-http/blob/v10.0.9/akka-http-core/src/main/scala/akka/http/impl/util/Rendering.scala#L206
       */
      val quotedString = "\"" + warning.replaceAll("\"","\\\\\"").replaceAll("[\\r\\n]+", " ").trim + "\""

      // https://www.w3.org/Protocols/rfc2616/rfc2616-sec14.html#sec14.46
      RawHeader("Warning", s"299 cromwell/$cromwellVersion $quotedString")
    }
  }
}

object WebServiceUtils extends WebServiceUtils {

  // Q: Why are we responding with a String of pretty-printed JSON instead of a JsObject?
  // A: There are customers who rely on the pretty printing to display the error directly in a terminal or GUI.
  // AEN 2018-12-05
  implicit class EnhancedThrowable(val e: Throwable) extends AnyVal {
    def failRequest(statusCode: StatusCode, warnings: Seq[String] = Vector.empty): Route = {
      completeResponse(statusCode, ContentTypes.`application/json`, prettyPrint(APIResponse.fail(e)), warnings)
    }
    def errorRequest(statusCode: StatusCode, warnings: Seq[String] = Vector.empty): Route = {
      completeResponse(statusCode, ContentTypes.`application/json`, prettyPrint(APIResponse.error(e)), warnings)
    }
  }

  private def prettyPrint(failureResponse: FailureResponse): String ={
    // .asJson cannot live inside a value class like `EnhancedThrowable`, hence the object method
    failureResponse.asJson.pretty(Printer.spaces2.copy(dropNullValues = true, colonLeft = ""))
  }

}
