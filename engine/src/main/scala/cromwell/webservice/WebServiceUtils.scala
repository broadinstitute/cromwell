package cromwell.webservice

import akka.http.scaladsl.model.{Multipart, StatusCode}
import akka.http.scaladsl.server.Route
import akka.stream.Materializer
import akka.util.{ByteString, Timeout}
import akka.http.scaladsl.marshalling.ToEntityMarshaller
import akka.http.scaladsl.model.headers.RawHeader
import akka.http.scaladsl.server.Directives.complete
import cromwell.webservice.routes.CromwellApiService._
import spray.json._
import cromwell.webservice.WorkflowJsonSupport.errorResponse

import scala.concurrent.{ExecutionContext, Future}

trait WebServiceUtils {

  type MaterializedFormData = Map[String, ByteString]

  def materializeFormData(formData: Multipart.FormData)(implicit timeout: Timeout, materializer: Materializer, executionContext: ExecutionContext): Future[MaterializedFormData] = {
    formData.parts.mapAsync[(String, ByteString)](1) {
      bodyPart => bodyPart.toStrict(timeout.duration)(materializer).map(strict => bodyPart.name -> strict.entity.data)(executionContext)
    }.runFold(Map.empty[String, ByteString])((map, tuple) => map + tuple)(materializer)
  }

  def completeResponse[A](statusCode: StatusCode, value: A, warnings: Seq[String])
                         (implicit mt: ToEntityMarshaller[A]): Route = {
    val warningHeaders = warnings.toIndexedSeq map { warning =>
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

    complete((statusCode, warningHeaders, value))
  }

}

object WebServiceUtils extends WebServiceUtils {

  implicit class EnhancedThrowable(val e: Throwable) extends AnyVal {
    def failRequest(statusCode: StatusCode, warnings: Seq[String] = Vector.empty): Route = {
      completeResponse(statusCode, APIResponse.fail(e).toJson.prettyPrint, warnings)
    }
    def errorRequest(statusCode: StatusCode, warnings: Seq[String] = Vector.empty): Route = {
      completeResponse(statusCode, APIResponse.error(e).toJson.prettyPrint, warnings)
    }
  }

}
