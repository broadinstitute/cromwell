package cromwell.engine.workflow.lifecycle.execution.callcaching

import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffActor.{CallDetails, HashDifference, SuccessfulCallCacheDiffResponse}
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import org.apache.commons.lang3.NotImplementedException
import spray.json._

object CallCacheDiffActorJsonFormatting extends SprayJsonSupport with DefaultJsonProtocol {

  implicit val callDetailsJsonFormatter = jsonFormat5(CallDetails)

  // Note: This json format is written out longform to get the non-standard Option behavior (the default omits 'None' fields altogether)
  implicit val hashDifferenceJsonFormatter = new RootJsonFormat[HashDifference] {
    override def write(hashDifference: HashDifference): JsValue = {
      def fromOption(opt: Option[String]) = opt.map(JsString.apply).getOrElse(JsNull)
      JsObject(Map(
        "hashKey" -> JsString(hashDifference.hashKey),
        "callA" -> fromOption(hashDifference.callA),
        "callB" -> fromOption(hashDifference.callB)
      ))
    }
    override def read(json: JsValue): HashDifference =
      throw new NotImplementedException("Programmer Error: No reader for HashDifferentials written. It was not expected to be required")
  }

  implicit val successfulResponseJsonFormatter = jsonFormat3(SuccessfulCallCacheDiffResponse)
}
