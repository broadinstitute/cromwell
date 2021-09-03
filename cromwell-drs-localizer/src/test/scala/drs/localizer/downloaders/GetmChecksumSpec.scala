package drs.localizer.downloaders

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks._
import cloud.nio.impl.drs.AccessUrl

class GetmChecksumSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  private val results = Table(
    ("Hashes", "Access URL", "expected"),
    (Option(Map("md5" -> "abcdefg")), "https://whatever", Md5("abcdefg")),
    (Option(Map("crc32c" -> "012345", "md5" -> "abcdefg")), "https://whatever", Md5("abcdefg")),
    (Option(Map("something weird" -> "012345", "md5" -> "abcdefg")), "https://whatever", Md5("abcdefg")),
    (Option(Map("something weird" -> "abcdefg", "crc32c" -> "012345")), "https://whatever", Crc32c("012345")),
    (Option(Map("etag" -> "abcdefg", "crc32c" -> "012345")), "https://whatever", Crc32c("012345")),
    (Option(Map("etag" -> "abcdefg", "something weird" -> "012345")), "https://whatever", Unsupported("etag, something weird")),
    (Option(Map("etag" -> "abcdefg", "something weird" -> "012345")), "https://whatever.s3.amazonaws.com/foo", AwsEtag("abcdefg")),
    (None, "https://whatever.s3.amazonaws.com/foo", Null),
  )

  forAll(results) { (hashes, url, expected) =>
    GetmChecksum(hashes, AccessUrl(url, None)) shouldBe expected
  }
}
