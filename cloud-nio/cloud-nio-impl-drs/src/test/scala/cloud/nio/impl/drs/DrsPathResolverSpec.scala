package cloud.nio.impl.drs

import java.time.{LocalDateTime, OffsetDateTime}

import cloud.nio.impl.drs.MarthaResponseSupport.convertMarthaResponseV2ToV3
import io.circe.{Json, JsonObject}
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

class DrsPathResolverSpec extends AnyFlatSpecLike with Matchers {
  private val mockGSA = SADataObject(data = Json.fromJsonObject(JsonObject("key"-> Json.fromString("value"))))
  private val crcHashValue = "8a366443"
  private val md5HashValue = "336ea55913bc261b72875bd259753046"
  private val shaHashValue = "f76877f8e86ec3932fd2ae04239fbabb8c90199dab0019ae55fa42b31c314c44"
  private val fullMarthaV2Response = MarthaV2Response(
    dos = DrsObject(
      data_object = DrsDataObject(
        size = Option(34905345),
        checksums = Option(Array(ChecksumObject(checksum = md5HashValue, `type` = "md5"), ChecksumObject(checksum = crcHashValue, `type` = "crc32c"))),
        updated = Option("2020-04-27T15:56:09.696Z"),
        urls = Array(Url("s3://my-s3-bucket/file-name"), Url("gs://my-gs-bucket/file-name"))
      )
    ),
    googleServiceAccount = Option(mockGSA)
  )

  private val fullMarthaV2ResponseNoTz = MarthaV2Response(
    dos = DrsObject(
      data_object = DrsDataObject(
        size = Option(34905345),
        checksums = Option(Array(ChecksumObject(checksum = md5HashValue, `type` = "md5"), ChecksumObject(checksum = crcHashValue, `type` = "crc32c"))),
        updated = Option("2020-01-15T17:46:25.694148"),
        urls = Array(Url("s3://my-s3-bucket/file-name"), Url("gs://my-gs-bucket/file-name"))
      )
    ),
    googleServiceAccount = Option(mockGSA)
  )

  private val fullMarthaResponse =  MarthaResponse(
    size = Option(34905345),
    timeUpdated = Option(OffsetDateTime.parse("2020-04-27T15:56:09.696Z").toString),
    bucket = Option("my-gs-bucket"),
    name = Option("file-name"),
    gsUri = Option("gs://my-gs-bucket/file-name"),
    googleServiceAccount = Option(mockGSA),
    hashes = Option(Map("md5" -> md5HashValue, "crc32c" -> crcHashValue))
  )

  private val fullMarthaResponseNoTz = MarthaResponse(
    size = Option(34905345),
    timeUpdated = Option(LocalDateTime.parse("2020-01-15T17:46:25.694148").toString),
    bucket = Option("my-gs-bucket"),
    name = Option("file-name"),
    gsUri = Option("gs://my-gs-bucket/file-name"),
    googleServiceAccount = Option(mockGSA),
    hashes = Option(Map("md5" -> md5HashValue, "crc32c" -> crcHashValue))
  )

  private val alfredHashValue = "xrd"
  private val completeHashesMap = Map(
    "betty" -> "abc123",
    "charles" -> "456",
    "alfred" -> alfredHashValue,
    "sha256" -> shaHashValue,
    "crc32c" -> crcHashValue,
    "md5" -> md5HashValue,
  )

  private val missingCRCHashesMap = Map(
    "alfred" -> alfredHashValue,
    "sha256" -> shaHashValue,
    "betty" -> "abc123",
    "md5" -> md5HashValue,
    "charles" -> "456",
  )

  private val onlySHAHashesMap = Map(
    "betty" -> "abc123",
    "charles" -> "456",
    "alfred" -> alfredHashValue,
    "sha256" -> shaHashValue,
  )

  private val noPreferredHashesMap = Map(
    "alfred" -> alfredHashValue,
    "betty" -> "abc123",
    "charles" -> "456",
  )

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

  behavior of "convertMarthaResponseV2ToV3()"

  it should "convert a full martha_v2 response to a the standard Martha response" in {
    convertMarthaResponseV2ToV3(fullMarthaV2Response) shouldBe fullMarthaResponse
  }

  it should "convert a full martha_v2 response to a the standard Martha response even if there is no timezone in `updated` field" in {
    convertMarthaResponseV2ToV3(fullMarthaV2ResponseNoTz) shouldBe fullMarthaResponseNoTz
  }
}
