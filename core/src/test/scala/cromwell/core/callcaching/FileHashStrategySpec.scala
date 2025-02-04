package cromwell.core.callcaching

import cromwell.core.callcaching.HashType.HashType
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

class FileHashStrategySpec extends AnyFlatSpecLike with Matchers {

  behavior of "FileHashStrategy.getFileHash"

  it should "return a hash when the strategy has one and one is supported by input" in {
    val strategy = FileHashStrategy.Crc32c
    val hash = "abc123"
    val returnedHash = strategy.getFileHash("my nice input",
                                            (input: String, hashType: HashType) =>
                                              hashType match {
                                                case HashType.Crc32c => Some(hash)
                                                case _ => None
                                              }
    )
    returnedHash shouldBe Some(FileHash(HashType.Crc32c, hash))
  }

  it should "return a hash when the strategy has one and many are supported by input" in {
    val strategy = FileHashStrategy.Crc32c
    val returnedHash = strategy.getFileHash(
      "my nice input",
      (input: String, hashType: HashType) =>
        hashType match {
          case HashType.Crc32c => Some("123")
          case HashType.Md5 => Some("456")
          case HashType.Etag => Some("789")
          case HashType.Sha256 => Some("1011")
          case HashType.Identity => Some(input)
          case _ => None
        }
    )
    returnedHash shouldBe Some(FileHash(HashType.Crc32c, "123"))
  }

  it should "return a hash when the strategy has many and many are supported by input" in {
    val strategy = FileHashStrategy.Drs
    val returnedHash = strategy.getFileHash(
      "my nice input",
      (input: String, hashType: HashType) =>
        hashType match {
          case HashType.Crc32c => Some("123")
          case HashType.Md5 => Some("456")
          case HashType.Etag => Some("789")
          case HashType.Sha256 => Some("1011")
          case HashType.Identity => Some(input)
          case _ => None
        }
    )
    returnedHash shouldBe Some(FileHash(HashType.Crc32c, "123"))
  }

  it should "fallback to the correct hash when the first one in the strategy is unavailable" in {
    val strategy = FileHashStrategy(List(HashType.Etag, HashType.Md5))
    val returnedHash = strategy.getFileHash(
      "my nice input",
      (input: String, hashType: HashType) =>
        hashType match {
          case HashType.Crc32c => Some("123")
          case HashType.Md5 => Some("456")
          case HashType.Sha256 => Some("1011")
          case HashType.Identity => Some(input)
          case _ => None
        }
    )
    returnedHash shouldBe Some(FileHash(HashType.Md5, "456"))
  }

  it should "return no hash when none of the defined strategies are available" in {
    val strategy = FileHashStrategy(List(HashType.Etag, HashType.Md5))
    val returnedHash = strategy.getFileHash(
      "my nice input",
      (input: String, hashType: HashType) =>
        hashType match {
          case HashType.Crc32c => Some("123")
          case HashType.Sha256 => Some("1011")
          case HashType.Identity => Some(input)
          case _ => None
        }
    )
    returnedHash shouldBe None
  }

  it should "not blow up when the strategy is empty" in {
    val strategy = FileHashStrategy.Empty
    val returnedHash = strategy.getFileHash(
      "my nice input",
      (input: String, hashType: HashType) =>
        hashType match {
          case HashType.Crc32c => Some("123")
          case HashType.Sha256 => Some("1011")
          case HashType.Identity => Some(input)
          case _ => None
        }
    )
    returnedHash shouldBe None
  }

  behavior of "FileHashStrategy.of"

  it should "construct a FileHashStrategy from an empty list" in {
    FileHashStrategy.of(List.empty) shouldBe FileHashStrategy(List.empty)
  }

  it should "construct a FileHashStrategy from an invalid string" in {
    FileHashStrategy.of(List("abc")) shouldBe FileHashStrategy(List.empty)
  }

  it should "construct a FileHashStrategy from a single valid string" in {
    FileHashStrategy.of(List("md5")) shouldBe FileHashStrategy(List(HashType.Md5))
  }

  it should "construct a FileHashStrategy from several valid strings" in {
    FileHashStrategy.of(List("MD5", "crc32c", "sHa256")) shouldBe FileHashStrategy(
      List(HashType.Md5, HashType.Crc32c, HashType.Sha256)
    )
  }

  it should "construct a FileHashStrategy from a mix of valid and invalid strings" in {
    FileHashStrategy.of(List("abc", "crc32c", "sdfalsdiwe", "etag")) shouldBe FileHashStrategy(
      List(HashType.Crc32c, HashType.Etag)
    )
  }

  behavior of "HashType.calculateHash"

  it should "calculate a crc32c hash" in {
    HashType.Crc32c.calculateHash("hello", true) shouldBe "9a71bb4c"
  }

  it should "calculate a base64-encoded crc32c hash" in {
    HashType.Crc32c.calculateHash("hello") shouldBe "mnG7TA=="
  }

  it should "calculate an etag hash on short data" in {
    HashType.Etag.calculateHash("hello") shouldBe "5d41402abc4b2a76b9719d911017c592"
  }

  it should "calculate an etag hash on medium data" in {
    val eightMB = 8 * 1024 * 1024
    val value = LazyList.continually(".").take(eightMB).mkString
    HashType.Etag.calculateHash(value) shouldBe "f89801f68b5028d64e0238ffb5a1b8e0"
  }

  it should "calculate an etag hash on long data" in {
    val eightMB = 8 * 1024 * 1024
    val value = LazyList.continually(".").take(eightMB + 1).mkString
    HashType.Etag.calculateHash(value) shouldBe "8e224b463f4f5202c9621820f7690a01-2"
  }

  it should "calculate an md5 hash" in {
    HashType.Md5.calculateHash("hello") shouldBe "5d41402abc4b2a76b9719d911017c592"
  }

  it should "calculate a sha256 hash" in {
    HashType.Sha256.calculateHash("hello") shouldBe "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824"
  }

  behavior of "FileHash.match"

  it should "correctly report a crc32c hex match" in {
    val hash = FileHash(HashType.Crc32c, "9a71bb4c")
    hash.matches("hello") shouldBe ChecksumSuccess()
  }

  it should "correctly report a crc32c b64 match" in {
    val hash = FileHash(HashType.Crc32c, "mnG7TA==")
    hash.matches("hello") shouldBe ChecksumSuccess()
  }

  it should "correctly report a crc32c mismatch" in {
    val hash = FileHash(HashType.Crc32c, "abc123")
    hash.matches("hello") shouldBe ChecksumFailure("9a71bb4c or mnG7TA==")
  }

  it should "correctly report an etag match on short data" in {
    val hash = FileHash(HashType.Etag, "5d41402abc4b2a76b9719d911017c592")
    hash.matches("hello") shouldBe ChecksumSuccess()
  }

  it should "correctly report an etag match on medium data" in {
    val eightMB = 8 * 1024 * 1024
    val value = LazyList.continually(".").take(eightMB).mkString
    val hash = FileHash(HashType.Etag, "f89801f68b5028d64e0238ffb5a1b8e0")
    hash.matches(value) shouldBe ChecksumSuccess()
  }

  it should "correctly report an etag match on long data" in {
    val eightMB = 8 * 1024 * 1024
    val value = LazyList.continually(".").take(eightMB + 1).mkString
    val hash = FileHash(HashType.Etag, "8e224b463f4f5202c9621820f7690a01-2")
    hash.matches(value) shouldBe ChecksumSuccess()
  }

  it should "correctly report an md5 match" in {
    val hash = FileHash(HashType.Md5, "5d41402abc4b2a76b9719d911017c592")
    hash.matches("hello") shouldBe ChecksumSuccess()
  }

  it should "correctly report a sha256 match" in {
    val hash = FileHash(HashType.Sha256, "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824")
    hash.matches("hello") shouldBe ChecksumSuccess()
  }

  it should "correctly report a sha256 mismatch" in {
    val hash = FileHash(HashType.Sha256, "not a very good SHA 256 hash")
    hash.matches("hello") shouldBe ChecksumFailure(
      "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824"
    )
  }
}
