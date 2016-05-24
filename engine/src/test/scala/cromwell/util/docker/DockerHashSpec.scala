package cromwell.util.docker

import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}

class DockerHashSpec extends FlatSpec with Matchers {
  behavior of "DockerHash"

  it should "parse valid digests" in {
    val digests = Table(
      ("digest", "hashType", "hashLength"),
      ("12345678", "unknown", 8),
      (":12345678", "", 8),
      ("type:12345678", "type", 8),
      ("type:11112222333344445555666677778888", "type", 32),
      ("type:1111222233334444555566667777888899990000aaaabbbbccccddddeeeeffff", "type", 64))

    forAll(digests) { (digest, hashType, hashLength) =>
      val dockerHash = DockerHash.fromDigest(digest).get
      dockerHash.hashType should be(hashType)
      dockerHash.hashString should have length hashLength
    }
  }

  it should "not parse invalid hash strings" in {
    val invalidHashStrings = Table(
      ("hashString", "exceptionMessage"),
      ("", "unexpected hash length: 0"),
      ("1", "hexBinary needs to be even-length: 1, unexpected hash length: 1"),
      ("hashhash", "contains illegal character for hexBinary: hashhash"),
      ("123456", "unexpected hash length: 6"),
      ("1234567890", "unexpected hash length: 10"))

    forAll(invalidHashStrings) { (hashString, exceptionMessage) =>
      val exception = DockerHash.fromHash("type", hashString).failed.get
      exception shouldBe an[IllegalArgumentException]
      exception.getMessage should be(s"hashString '$hashString' is not valid: $exceptionMessage")
    }
  }

  it should "create hash digests" in {
    DockerHash("type", "12345678").digest should be("type:12345678")
  }

  it should "not rehash an empty hash" in {
    val emptyHash = Vector.empty[DockerHash]
    val exception = DockerHash.fromSeq("collection", emptyHash).failed.get
    exception shouldBe an[IllegalArgumentException]
    exception.getMessage should be("docker hashes is empty")
  }

  it should "rehash a sequence of hashes" in {
    val hash1 = DockerHash.fromHash("single", "12345678")
    val hash2 = DockerHash.fromHash("single", "87654321")
    val hashes = Vector(hash1, hash2)
    val rehash = DockerHash.fromTries("collection", hashes).get
    rehash.hashType should be("collection-single-md5")
    rehash.hashString should be("461595e4bdb090ce41e7818287954d86")
  }

  it should "not rehash a sequence of different hash types" in {
    val hash1 = DockerHash.fromHash("hash1", "12345678")
    val hash2 = DockerHash.fromHash("hash2", "12345678")
    val hashes = Vector(hash1, hash2)
    val exception = DockerHash.fromTries("collection", hashes).failed.get
    exception shouldBe an[IllegalArgumentException]
    exception.getMessage should be("found more than one docker hash type: Vector(hash1, hash2)")
  }
}
