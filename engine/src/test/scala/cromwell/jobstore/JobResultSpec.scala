package cromwell.jobstore

import cromwell.jobstore.JobResultJsonFormatter._
import org.scalatest.{FlatSpec, Matchers}
import spray.json._
import wom.JobOutput
import wom.types.{WdlIntegerType, WdlMapType, WdlStringType}
import wom.values.{WdlInteger, WdlMap, WdlString}

class JobResultSpec extends FlatSpec with Matchers {

  behavior of "JobResult Json Serialization/Deserialization"

  it should "write and read JSON for Job successes" in {

    val success = JobResultSuccess(Some(0), Map("abc" -> JobOutput(WdlString("hello"))))
    val asJson = success.toJson

    val jsonString = asJson.toString()
    jsonString shouldBe """{"returnCode":0,"jobOutputs":{"abc":"hello"}}"""

    val fromJsonString = jsonString.parseJson
    val fromJson = fromJsonString.convertTo[JobResultSuccess]

    fromJson shouldBe success
  }

  it should "write and read more complicated WdlValues" in {
    val success = JobResultSuccess(Some(0), Map("abc" -> JobOutput(WdlMap(WdlMapType(WdlStringType, WdlIntegerType), Map(WdlString("hello") -> WdlInteger(4), WdlString("goodbye") -> WdlInteger(6))))))
    val asJson = success.toJson

    val jsonString = asJson.toString()
    val fromJsonString = jsonString.parseJson
    val fromJson = fromJsonString.convertTo[JobResultSuccess]

    fromJson shouldBe success
  }

  it should "write and read JSON for job failures" in {

    val failure = JobResultFailure(Some(0), new Exception("abc"), retryable = false)
    val asJson = failure.toJson

    val jsonString = asJson.toString()
    jsonString shouldBe """{"returnCode":0,"reason":"abc","retryable":false}"""

    val fromJsonString = jsonString.parseJson
    val fromJson = fromJsonString.convertTo[JobResultFailure]

    fromJson.returnCode shouldBe failure.returnCode
    fromJson.reason.getMessage shouldBe failure.reason.getMessage
  }
}
