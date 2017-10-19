package cromwell.jobstore

import cromwell.jobstore.JobResultJsonFormatter._
import cromwell.util.WomMocks
import org.scalatest.{FlatSpec, Matchers}
import spray.json._
import wom.types.{WomIntegerType, WomMapType, WomStringType}
import wom.values.{WomInteger, WomMap, WomString}

class JobResultSpec extends FlatSpec with Matchers {

  behavior of "JobResult Json Serialization"

  it should "write JSON for Job successes" in {

    val success = JobResultSuccess(Some(0), WomMocks.mockOutputExpectations(Map("abc" -> WomString("hello"))))
    val asJson = success.toJson

    val jsonString = asJson.toString()
    jsonString shouldBe """{"returnCode":0,"jobOutputs":{"abc":"hello"}}"""
  }

  it should "write more complicated WdlValues" in {
    val success = JobResultSuccess(Some(0), WomMocks.mockOutputExpectations(Map("abc" -> WomMap(WomMapType(WomStringType, WomIntegerType), Map(WomString("hello") -> WomInteger(4), WomString("goodbye") -> WomInteger(6))))))
    val asJson = success.toJson

    val jsonString = asJson.toString()
    jsonString shouldBe """{"returnCode":0,"jobOutputs":{"abc":{"hello":4,"goodbye":6}}}"""
  }

  it should "write JSON for job failures" in {

    val failure = JobResultFailure(Some(0), new Exception("abc"), retryable = false)
    val asJson = failure.toJson

    val jsonString = asJson.toString()
    jsonString shouldBe """{"returnCode":0,"reason":"abc","retryable":false}"""
  }
}
