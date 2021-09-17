package drs.localizer.downloaders

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks._
import cloud.nio.impl.drs.AccessUrl

class GetmChecksumSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  it should "deduce the correct getm checksum algorithm" in {
    val results = Table(
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

  {
    val results = Table(
      ("description", "algorithm", "expected"),
      ("md5", Md5("abcdef"), "--checksum-algorithm 'md5' --checksum abcdef"),
      ("crc32c", Crc32c("012345"), "--checksum-algorithm 'gs_crc32c' --checksum 012345"),
      ("AWS ETag", AwsEtag("012345"), "--checksum-algorithm 's3_etag' --checksum 012345"),
      // Escape checksum values constructed from unvalidated data returned by DRS servers.
      ("Unsupported", Unsupported("Robert'); DROP TABLE Students; --"), raw"--checksum-algorithm 'null' --checksum Robert\'\)\;\ DROP\ TABLE\ Students\;\ --"),
      ("Null", Null, "--checksum-algorithm 'null' --checksum null"),
    )

    forAll(results) { (description, algorithm, expected) =>
      it should s"produce the expected checksum arguments for $description" in {
        algorithm.args shouldBe (expected)
      }
    }
  }
}
