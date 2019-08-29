package cromwell.engine.workflow.lifecycle.execution.callcaching

import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffActor.{CallDetails, HashDifference, SuccessfulCallCacheDiffResponse}
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import spray.json._

object CallCacheDiffActorJsonFormatting extends SprayJsonSupport with DefaultJsonProtocol {

  implicit val callDetailsJsonFormatter = jsonFormat5(CallDetails)
  implicit val hashDifferenceJsonFormatter = jsonFormat3(HashDifference)
  implicit val successfulResponseJsonFormatter = jsonFormat3(SuccessfulCallCacheDiffResponse)
}
