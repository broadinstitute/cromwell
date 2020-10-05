package cloud.nio.impl.drs

import java.nio.file.attribute.FileTime
import java.time.OffsetDateTime

import cloud.nio.impl.drs.DrsCloudNioRegularFileAttributes._
import cloud.nio.impl.drs.MarthaResponseSupport.convertMarthaResponseV2ToV3
import io.circe.{Json, JsonObject}
import org.apache.http.ProtocolVersion
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import shapeless.{Lens, lens}

class DrsPathResolverSpec extends AnyFlatSpecLike with Matchers {
  private val mockGSA = SADataObject(data = Json.fromJsonObject(JsonObject("key"-> Json.fromString("value"))))
  private val crcHashValue = "8a366443"
  private val md5HashValue = "336ea55913bc261b72875bd259753046"
  private val shaHashValue = "f76877f8e86ec3932fd2ae04239fbabb8c90199dab0019ae55fa42b31c314c44"

  private val marthaV2UpdatedLens =
    lens[MarthaV2Response]
      .dos
      .data_object
      .updated
      // https://stackoverflow.com/questions/23874998/shapeless-lenses-in-idea#23899094
      .asInstanceOf[Lens[MarthaV2Response, Option[String]]]

  private val fullMarthaV2Response = MarthaV2Response(
    dos = DosObject(
      data_object = DosDataObject(
        name = Option("actual_file_name"),
        size = Option(34905345),
        checksums = Option(Array(ChecksumObject(checksum = md5HashValue, `type` = "md5"), ChecksumObject(checksum = crcHashValue, `type` = "crc32c"))),
        created = Option("2020-04-27T15:56:09.696Z"),
        updated = Option("2020-04-27T15:56:09.696Z"),
        urls = Array(Url("s3://my-s3-bucket/file-name"), Url("gs://my-gs-bucket/file-name"))
      )
    ),
    googleServiceAccount = Option(mockGSA)
  )

  private val fullMarthaV2ResponseNoTz =
    marthaV2UpdatedLens
      .set(fullMarthaV2Response)(fullMarthaV2Response.dos.data_object.updated.map(_.stripSuffix("Z")))

  private val fullMarthaResponse = MarthaResponse(
    size = Option(34905345),
    timeCreated = Option(OffsetDateTime.parse("2020-04-27T15:56:09.696Z").toString),
    timeUpdated = Option(OffsetDateTime.parse("2020-04-27T15:56:09.696Z").toString),
    gsUri = Option("gs://my-gs-bucket/file-name"),
    googleServiceAccount = Option(mockGSA),
    fileName = Option("actual_file_name"),
    hashes = Option(Map("md5" -> md5HashValue, "crc32c" -> crcHashValue))
  )

  private val fullMarthaResponseNoTz =
    fullMarthaResponse
      .copy(timeUpdated = fullMarthaResponse.timeUpdated.map(_.stripSuffix("Z")))

  private val fullMarthaResponseNoTime =
    fullMarthaResponse
      .copy(timeUpdated = None)

  private val fullMarthaResponseBadTz =
    fullMarthaResponse
      .copy(timeUpdated = fullMarthaResponse.timeUpdated.map(_.stripSuffix("Z") + "BADTZ"))

  private val alfredHashValue = "xrd"
  private val completeHashesMap = Option(Map(
    "betty" -> "abc123",
    "charles" -> "456",
    "alfred" -> alfredHashValue,
    "sha256" -> shaHashValue,
    "crc32c" -> crcHashValue,
    "md5" -> md5HashValue,
  ))

  private val missingCRCHashesMap = Option(Map(
    "alfred" -> alfredHashValue,
    "sha256" -> shaHashValue,
    "betty" -> "abc123",
    "md5" -> md5HashValue,
    "charles" -> "456",
  ))

  private val onlySHAHashesMap = Option(Map(
    "betty" -> "abc123",
    "charles" -> "456",
    "alfred" -> alfredHashValue,
    "sha256" -> shaHashValue,
  ))

  private val noPreferredHashesMap = Option(Map(
    "alfred" -> alfredHashValue,
    "betty" -> "abc123",
    "charles" -> "456",
  ))

  behavior of "fileHash()"

  it should "return crc32c hash from `hashes` in Martha response when there is a crc32c" in {
    DrsCloudNioRegularFileAttributes.getPreferredHash(completeHashesMap) shouldBe Option(crcHashValue)
  }

  it should "return md5 hash from `hashes` in Martha response when there is no crc32c" in {
    DrsCloudNioRegularFileAttributes.getPreferredHash(missingCRCHashesMap) shouldBe Option(md5HashValue)
  }

  it should "return sha256 hash from `hashes` in Martha response when there is only a sha256" in {
    DrsCloudNioRegularFileAttributes.getPreferredHash(onlySHAHashesMap) shouldBe Option(shaHashValue)
  }

  it should "return first (alphabetized by type) hash from `hashes` in Martha response when there are no preferred hash types" in {
    DrsCloudNioRegularFileAttributes.getPreferredHash(noPreferredHashesMap) shouldBe Option(alfredHashValue)
  }

  it should "return None when no hashes object is returned" in {
    DrsCloudNioRegularFileAttributes.getPreferredHash(None) shouldBe Option(None)
  }

  it should "return None when an empty hash object is returned" in {
    DrsCloudNioRegularFileAttributes.getPreferredHash(Option(Map.empty)) shouldBe Option(None)
  }

  behavior of "convertMarthaResponseV2ToV3()"

  it should "convert a full martha_v2 response to a the standard Martha response" in {
    convertMarthaResponseV2ToV3(fullMarthaV2Response) shouldBe fullMarthaResponse
  }

  it should "convert a full martha_v2 response to a the standard Martha response even if there is no timezone in `updated` field" in {
    convertMarthaResponseV2ToV3(fullMarthaV2ResponseNoTz) shouldBe fullMarthaResponseNoTz
  }

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
    import cloud.nio.impl.drs.MarthaResponseSupport.marthaFailureResponseDecoder

    val maybeDecoded = decode[MarthaFailureResponse](failureResponseJson)
    maybeDecoded map { decoded: MarthaFailureResponse =>
      decoded.response.text shouldBe "{\"msg\":\"User 'null' does not have required action: read_data\",\"status_code\":500}"
    }
  }

  import org.apache.http.message.BasicStatusLine

  val drsPathForDebugging = "drs://my_awesome_drs"
  val responseStatusLine = new BasicStatusLine(new ProtocolVersion("http", 1, 2) , 345, "test-reason")
  val testMarthaUri = "www.martha_v3.com"

  it should "construct an error message from a populated, well-formed failure response" in {
    val failureResponse = Option(failureResponseJson)

    MarthaResponseSupport.errorMessageFromResponse(drsPathForDebugging, failureResponse, responseStatusLine, testMarthaUri) shouldBe {
      "Could not access object 'drs://my_awesome_drs'. Status: 345, reason: 'test-reason', Martha location: 'www.martha_v3.com', message: '{\"msg\":\"User 'null' does not have required action: read_data\",\"status_code\":500}'"
    }
  }

  it should "construct an error message from an empty failure response" in {
    MarthaResponseSupport.errorMessageFromResponse(drsPathForDebugging, None, responseStatusLine, testMarthaUri) shouldBe {
      "Could not access object 'drs://my_awesome_drs'. Status: 345, reason: 'test-reason', Martha location: 'www.martha_v3.com', message: (empty response)"
    }
  }

  // Technically we enter this case when preparing the "error message" for a successful response, because naturally `MarthaResponse` does not deserialize to `MarthaFailureResponse`
  // But then there's no error so we throw it away :shrug:
  it should "construct an error message from a malformed failure response" in {
    val unparsableFailureResponse = Option("something went horribly wrong")

    MarthaResponseSupport.errorMessageFromResponse(drsPathForDebugging, unparsableFailureResponse, responseStatusLine, testMarthaUri) shouldBe {
      "Could not access object 'drs://my_awesome_drs'. Status: 345, reason: 'test-reason', Martha location: 'www.martha_v3.com', message: 'something went horribly wrong'"
    }
  }

  it should "resolve an ISO-8601 date with timezone" in {
    val lastModifiedTimeIO = convertToFileTime(
      "drs://my_awesome_drs",
      MarthaField.TimeUpdated,
      fullMarthaResponse.timeUpdated,
    )
    lastModifiedTimeIO.unsafeRunSync() should
      be(Option(FileTime.from(OffsetDateTime.parse("2020-04-27T15:56:09.696Z").toInstant)))
  }

  it should "resolve an ISO-8601 date without timezone" in {
    val lastModifiedTimeIO = convertToFileTime(
      "drs://my_awesome_drs",
      MarthaField.TimeUpdated,
      fullMarthaResponseNoTz.timeUpdated,
    )
    lastModifiedTimeIO.unsafeRunSync() should
      be(Option(FileTime.from(OffsetDateTime.parse("2020-04-27T15:56:09.696Z").toInstant)))
  }

  it should "not resolve an date that does not contain a timeUpdated" in {
    val lastModifiedTimeIO = convertToFileTime(
      "drs://my_awesome_drs",
      MarthaField.TimeUpdated,
      fullMarthaResponseNoTime.timeUpdated,
    )
    lastModifiedTimeIO.unsafeRunSync() should be(None)
  }

  it should "not resolve an date that is not ISO-8601" in {
    val lastModifiedTimeIO = convertToFileTime(
      "drs://my_awesome_drs",
      MarthaField.TimeUpdated,
      fullMarthaResponseBadTz.timeUpdated,
    )
    the[RuntimeException] thrownBy lastModifiedTimeIO.unsafeRunSync() should have message
      "Error while parsing 'timeUpdated' value from Martha to FileTime for DRS path drs://my_awesome_drs. " +
        "Reason: DateTimeParseException: Text '2020-04-27T15:56:09.696BADTZ' could not be parsed at index 23."
  }
}
