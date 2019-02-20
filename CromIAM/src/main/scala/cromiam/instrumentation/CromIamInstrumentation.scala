package cromiam.instrumentation

import akka.http.scaladsl.model.{HttpRequest, HttpResponse}
import cats.data.NonEmptyList
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

  val samPrefix: NonEmptyList[String] = NonEmptyList.one("sam")
  val getWhitelistPrefix = NonEmptyList.one("get-whitelist")
  val userCollectionPrefix = NonEmptyList.one("user-collection")
  val authCollectionPrefix = NonEmptyList.one("auth-collection")
  val registerCollectionPrefix = NonEmptyList.one("register-collection")
  val rootWfIdPrefix = NonEmptyList.one("root-workflow-id")
  val wfCollectionPrefix = NonEmptyList.one("workflow-collection")


  def convertRequestToPath(httpRequest: HttpRequest): NonEmptyList[String] = NonEmptyList.of(
    // Returns the path of the URI only, without query parameters (e.g: api/engine/workflows/metadata)
    httpRequest.uri.path.toString().stripPrefix("/")
      // Replace UUIDs with [id] to keep paths same regardless of the workflow
      .replaceAll(CromIamInstrumentation.UUIDRegex, "[id]"),
    // Name of the method (e.g: GET)
    httpRequest.method.value
  )

  def makePathFromRequestAndResponseString(httpRequest: HttpRequest, response: String): InstrumentationPath =
    convertRequestToPath(httpRequest).concatNel(NonEmptyList.of(response))

  def makePathFromRequestAndResponse(httpRequest: HttpRequest, httpResponse: HttpResponse): InstrumentationPath =
    convertRequestToPath(httpRequest).concatNel(NonEmptyList.of(httpResponse.status.intValue.toString))

  def sendTimingApi(statsDPath: InstrumentationPath, timing: FiniteDuration, prefixToStatsd: NonEmptyList[String]): Unit = {
    sendTiming(prefixToStatsd.concatNel(statsDPath), timing, CromIamPrefix)
  }

  def instrumentationPrefixForSam(methodPrefix: NonEmptyList[String]): NonEmptyList[String] = samPrefix.concatNel(methodPrefix)

  def instrumentRequest[A](func: () => FailureResponseOrT[A],
                           httpRequest: HttpRequest,
                           prefix: NonEmptyList[String]): FailureResponseOrT[A] = {
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
