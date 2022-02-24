package cromiam.instrumentation

import akka.http.scaladsl.model.{HttpRequest, HttpResponse}
import cromwell.api.model.FailureResponseOrT
import cromwell.core.instrumentation.InstrumentationPrefixes._
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

  val samPrefix: InstrumentationPath = InstrumentationPath.withParts("sam")
  val getWhitelistPrefix: InstrumentationPath = InstrumentationPath.withParts("get-whitelist")
  val userCollectionPrefix: InstrumentationPath = InstrumentationPath.withParts("user-collection")
  val authCollectionPrefix: InstrumentationPath = InstrumentationPath.withParts("auth-collection")
  val registerCollectionPrefix: InstrumentationPath = InstrumentationPath.withParts("register-collection")
  val rootWfIdPrefix: InstrumentationPath = InstrumentationPath.withParts("root-workflow-id")
  val wfCollectionPrefix: InstrumentationPath = InstrumentationPath.withParts("workflow-collection")


  def convertRequestToPath(httpRequest: HttpRequest): InstrumentationPath = InstrumentationPath
    // Returns the path of the URI only, without query parameters (e.g: api/engine/workflows/metadata)
    .withHighVariantPart("uri", httpRequest.uri.path.toString().stripPrefix("/")
      // Replace UUIDs with [id] to keep paths same regardless of the workflow
      .replaceAll(CromIamInstrumentation.UUIDRegex, "[id]"))
    // Name of the method (e.g: GET)
    .withHighVariantPart("method", httpRequest.method.value)

  def makePathFromRequestAndResponseString(httpRequest: HttpRequest, response: String): InstrumentationPath =
    convertRequestToPath(httpRequest).withHighVariantPart("response", response)

  def makePathFromRequestAndResponse(httpRequest: HttpRequest, httpResponse: HttpResponse): InstrumentationPath =
    convertRequestToPath(httpRequest).withHighVariantPart("response", httpResponse.status.intValue().toString)

  def sendTimingApi(statsDPath: InstrumentationPath, timing: FiniteDuration, prefixToStatsd: InstrumentationPath): Unit = {
    sendTiming(prefixToStatsd.concat(statsDPath), timing, CromIamPrefix)
  }

  def instrumentationPrefixForSam(methodPrefix: InstrumentationPath): InstrumentationPath = samPrefix.concat(methodPrefix)

  def instrumentRequest[A](func: () => FailureResponseOrT[A],
                           httpRequest: HttpRequest,
                           prefix: InstrumentationPath): FailureResponseOrT[A] = {
    def now(): Deadline = Deadline.now

    val startTimestamp = now()

    def sendDuration(instrumentationPath: InstrumentationPath): Unit = {
      val duration = now() - startTimestamp
      sendTimingApi(instrumentationPath, duration, prefix)
    }

    func() map {
      case response: HttpResponse =>
        sendDuration(makePathFromRequestAndResponse(httpRequest, response))
        response.asInstanceOf[A]
      case other =>
        sendDuration(makePathFromRequestAndResponseString(httpRequest, "success"))
        other.asInstanceOf[A]
    } leftMap { error =>
      sendDuration(makePathFromRequestAndResponseString(httpRequest, "failure"))
      error
    }
  }
}
