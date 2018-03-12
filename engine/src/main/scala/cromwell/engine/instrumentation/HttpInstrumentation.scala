package cromwell.engine.instrumentation

import akka.http.scaladsl.model.{HttpRequest, HttpResponse}
import akka.http.scaladsl.server.Directive0
import akka.http.scaladsl.server.Directives.{extractRequest, mapResponse}
import cats.data.NonEmptyList
import cromwell.core.instrumentation.InstrumentationPrefixes._
import cromwell.services.instrumentation.{CromwellInstrumentation, CromwellInstrumentationActor}
import cromwell.services.instrumentation.CromwellInstrumentation._

import scala.concurrent.duration._

object HttpInstrumentation {
  private val UUIDRegex = {
    val hex = """[a-fA-F\d]"""
    s"""$hex{8}-$hex{4}-$hex{4}-$hex{4}-$hex{12}"""
  }
}

trait HttpInstrumentation extends CromwellInstrumentation {
  
  private def makeRequestPath(httpRequest: HttpRequest, httpResponse: HttpResponse): InstrumentationPath = NonEmptyList.of(
    // Returns the path of the URI only, without query parameters (e.g: api/engine/workflows/metadata)
    httpRequest.uri.path.toString().stripPrefix("/")
      // Replace UUIDs with [id] to keep paths unique regardless of the workflow
      .replaceAll(HttpInstrumentation.UUIDRegex, "[id]"),
    // Name of the method (e.g: GET)
    httpRequest.method.value,
    // Status code of the Response (e.g: 200)
    httpResponse.status.intValue.toString
  )

  private def sendTimingApi(statsDPath: InstrumentationPath, timing: FiniteDuration) = {
    sendTiming(statsDPath, timing, ApiPrefix)
  }

  def instrumentRequest: Directive0 = extractRequest flatMap { request =>
    val timeStamp = System.currentTimeMillis
    mapResponse { response =>
      /*
        * Send a metric corresponding to the request response time.
        * Note: The current StatsD implementation always pairs a timing metric with a counter metric
        * So no need to explicitly send one
       */
      sendTimingApi(makeRequestPath(request, response), (System.currentTimeMillis - timeStamp).millis)
      response
    }
  }
}
