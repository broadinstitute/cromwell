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
}
