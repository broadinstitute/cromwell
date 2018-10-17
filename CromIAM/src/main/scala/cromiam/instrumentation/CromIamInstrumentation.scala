package cromiam.instrumentation

import akka.http.scaladsl.model.HttpRequest
import cats.data.NonEmptyList
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.instrumentation.CromwellInstrumentation.InstrumentationPath

import scala.concurrent.duration._

object CromIamInstrumentation {
  private val UUIDRegex = {
    val hex = """[a-fA-F\d]"""
    s"""$hex{8}-$hex{4}-$hex{4}-$hex{4}-$hex{12}"""
  }
}

trait CromIamInstrumentation extends CromwellInstrumentation {

  def makeRequestPath(httpRequest: HttpRequest, response: String): InstrumentationPath = NonEmptyList.of(
    // Returns the path of the URI only, without query parameters (e.g: api/engine/workflows/metadata)
    httpRequest.uri.path.toString().stripPrefix("/")
      // Replace UUIDs with [id] to keep paths unique regardless of the workflow
      .replaceAll(CromIamInstrumentation.UUIDRegex, "[id]"),
    // Name of the method (e.g: GET)
    httpRequest.method.value,
    // Status code of the Response (e.g: 200)
    response
  )

  def sendTimingApi(statsDPath: InstrumentationPath, timing: FiniteDuration, prefix: String) = {
    val prefixB = Option(s"cromiam.sam.$prefix")
    sendTiming(statsDPath, timing, prefixB)
  }

//  def instrumentRequest: Directive0 = extractRequest flatMap { request =>
//    val timeStamp = System.currentTimeMillis
//    mapResponse { response =>
//      /*
//        * Send a metric corresponding to the request response time.
//        * Note: The current StatsD implementation always pairs a timing metric with a counter metric
//        * So no need to explicitly send one
//       */
//      sendTimingApi(makeRequestPath(request, response), (System.currentTimeMillis - timeStamp).millis)
//      response
//    }
//  }
}
