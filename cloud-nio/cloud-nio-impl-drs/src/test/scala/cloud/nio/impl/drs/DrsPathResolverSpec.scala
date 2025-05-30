package cloud.nio.impl.drs

import cats.data.NonEmptyList

import java.nio.file.attribute.FileTime
import java.time.OffsetDateTime
import cloud.nio.impl.drs.DrsCloudNioRegularFileAttributes._
import common.assertion.CromwellTimeoutSpec
import io.circe.{Json, JsonObject}
import org.apache.http.ProtocolVersion
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

class DrsPathResolverSpec extends AnyFlatSpecLike with CromwellTimeoutSpec with Matchers {
  private val mockGSA = SADataObject(data = Json.fromJsonObject(JsonObject("key" -> Json.fromString("value"))))
  private val crcHashValue = "8a366443"
  private val md5HashValue = "336ea55913bc261b72875bd259753046"

  private val fullDrsResolverResponse = DrsResolverResponse(
    size = Option(34905345),
    timeCreated = Option(OffsetDateTime.parse("2020-04-27T15:56:09.696Z").toString),
    timeUpdated = Option(OffsetDateTime.parse("2020-04-27T15:56:09.696Z").toString),
    gsUri = Option("gs://my-gs-bucket/file-name"),
    googleServiceAccount = Option(mockGSA),
    fileName = Option("actual_file_name"),
    hashes = Option(Map("md5" -> md5HashValue, "crc32c" -> crcHashValue))
  )

  private val fullDrsResolverResponseNoTz =
    fullDrsResolverResponse
      .copy(timeUpdated = fullDrsResolverResponse.timeUpdated.map(_.stripSuffix("Z")))

  private val fullDrsResolverResponseNoTime =
    fullDrsResolverResponse
      .copy(timeUpdated = None)

  private val fullDrsResolverResponseBadTz =
    fullDrsResolverResponse
      .copy(timeUpdated = fullDrsResolverResponse.timeUpdated.map(_.stripSuffix("Z") + "BADTZ"))

  private val failureResponseJson = """
    {
      "status": 500,
      "response": {
        "req": {
        "method": "GET",
        "url": "https://jade.datarepo-dev.broadinstitute.org/ga4gh/drs/v1/objects/v1_0c86170e-312d-4b39-a0a4-2a2bfaa24c7a_c0e40912-8b14-43f6-9a2f-b278144d0060",
        "headers": {
        "user-agent": "node-superagent/3.8.3",
        "authori - not gonna get me this time, git-secrets": "A bear with a token"
      }
      },
        "header": {
        "date": "Wed, 09 Sep 2020 14:52:10 GMT",
        "server": "nginx/1.18.0",
        "x-frame-options": "SAMEORIGIN",
        "content-type": "application/json;charset=UTF-8",
        "transfer-encoding": "chunked",
        "via": "1.1 google",
        "alt-svc": "clear",
        "connection": "close"
      },
        "status": 500,
        "text": "{\"msg\":\"User 'null' does not have required action: read_data\",\"status_code\":500}"
      }
    }"""

  it should "successfully parse a failure" in {
    import io.circe.parser.decode
    import cloud.nio.impl.drs.DrsResolverResponseSupport.drsResolverFailureResponseDecoder

    val maybeDecoded = decode[DrsResolverFailureResponse](failureResponseJson)
    maybeDecoded map { decoded: DrsResolverFailureResponse =>
      decoded.response.text shouldBe "{\"msg\":\"User 'null' does not have required action: read_data\",\"status_code\":500}"
    }
  }

  import org.apache.http.message.BasicStatusLine

  val drsPathForDebugging = "drs://my_awesome_drs"
  val responseStatusLine = new BasicStatusLine(new ProtocolVersion("http", 1, 2), 345, "test-reason")
  val testDrsResolverUri = "www.drshub_v4.com"

  it should "construct the right request when using Azure creds" in {
    val resolver = new MockDrsPathResolver(drsCredentials = AzureDrsCredentials())
    val drsRequest = resolver.makeDrsResolverRequest(drsPathForDebugging, NonEmptyList.of(DrsResolverField.AccessUrl))
    drsRequest.cloudPlatform shouldBe Option(DrsCloudPlatform.Azure)
  }

  it should "construct the right request when using Google creds" in {
    val resolver = new MockDrsPathResolver()
    val drsRequest = resolver.makeDrsResolverRequest(drsPathForDebugging, NonEmptyList.of(DrsResolverField.AccessUrl))
    drsRequest.cloudPlatform shouldBe Option(DrsCloudPlatform.GoogleStorage)
  }

  it should "construct an error message from a populated, well-formed failure response" in {
    val failureResponse = Option(failureResponseJson)

    DrsResolverResponseSupport.errorMessageFromResponse(drsPathForDebugging,
                                                        failureResponse,
                                                        responseStatusLine,
                                                        testDrsResolverUri
    ) shouldBe {
      "Could not access object 'drs://my_awesome_drs'. Status: 345, reason: 'test-reason', DRS Resolver location: 'www.drshub_v4.com', message: '{\"msg\":\"User 'null' does not have required action: read_data\",\"status_code\":500}'"
    }
  }

  it should "construct an error message from an empty failure response" in {
    DrsResolverResponseSupport.errorMessageFromResponse(drsPathForDebugging,
                                                        None,
                                                        responseStatusLine,
                                                        testDrsResolverUri
    ) shouldBe {
      "Could not access object 'drs://my_awesome_drs'. Status: 345, reason: 'test-reason', DRS Resolver location: 'www.drshub_v4.com', message: (empty response)"
    }
  }

  // Technically we enter this case when preparing the "error message" for a successful response, because naturally `DrsResolverResponse` does not deserialize to `DrsResolverFailureResponse`
  // But then there's no error so we throw it away :shrug:
  it should "construct an error message from a malformed failure response" in {
    val unparsableFailureResponse = Option("something went horribly wrong")

    DrsResolverResponseSupport.errorMessageFromResponse(drsPathForDebugging,
                                                        unparsableFailureResponse,
                                                        responseStatusLine,
                                                        testDrsResolverUri
    ) shouldBe {
      "Could not access object 'drs://my_awesome_drs'. Status: 345, reason: 'test-reason', DRS Resolver location: 'www.drshub_v4.com', message: 'something went horribly wrong'"
    }
  }

  it should "resolve an ISO-8601 date with timezone" in {
    val lastModifiedTimeIO = convertToFileTime(
      "drs://my_awesome_drs",
      DrsResolverField.TimeUpdated,
      fullDrsResolverResponse.timeUpdated
    )
    lastModifiedTimeIO.unsafeRunSync() should
      be(Option(FileTime.from(OffsetDateTime.parse("2020-04-27T15:56:09.696Z").toInstant)))
  }

  it should "resolve an ISO-8601 date without timezone" in {
    val lastModifiedTimeIO = convertToFileTime(
      "drs://my_awesome_drs",
      DrsResolverField.TimeUpdated,
      fullDrsResolverResponseNoTz.timeUpdated
    )
    lastModifiedTimeIO.unsafeRunSync() should
      be(Option(FileTime.from(OffsetDateTime.parse("2020-04-27T15:56:09.696Z").toInstant)))
  }

  it should "not resolve an date that does not contain a timeUpdated" in {
    val lastModifiedTimeIO = convertToFileTime(
      "drs://my_awesome_drs",
      DrsResolverField.TimeUpdated,
      fullDrsResolverResponseNoTime.timeUpdated
    )
    lastModifiedTimeIO.unsafeRunSync() should be(None)
  }

  it should "not resolve an date that is not ISO-8601" in {
    val lastModifiedTimeIO = convertToFileTime(
      "drs://my_awesome_drs",
      DrsResolverField.TimeUpdated,
      fullDrsResolverResponseBadTz.timeUpdated
    )
    the[RuntimeException] thrownBy lastModifiedTimeIO.unsafeRunSync() should have message
      "Error while parsing 'timeUpdated' value from DRS Resolver to FileTime for DRS path drs://my_awesome_drs. " +
      "Reason: DateTimeParseException: Text '2020-04-27T15:56:09.696BADTZ' could not be parsed at index 23."
  }
}
