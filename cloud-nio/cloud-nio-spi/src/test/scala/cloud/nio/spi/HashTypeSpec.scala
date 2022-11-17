package cloud.nio.spi

import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

class HashTypeSpec extends AnyFlatSpecLike with Matchers {

  behavior of "HashType.calculateHash"

  it should "calculate a crc32c hash" in {
    HashType.Crc32c.calculateHash("hello") shouldBe "9a71bb4c"
  }

  it should "calculate a GCS crc32c hash" in {
    HashType.GcsCrc32c.calculateHash("hello") shouldBe "mnG7TA=="
  }

  it should "calculate an etag hash on short data" in {
    HashType.S3Etag.calculateHash("hello") shouldBe "5d41402abc4b2a76b9719d911017c592"
  }

  it should "calculate an etag hash on medium data" in {
    val eightMB = 8 * 1024 * 1024
    val value = LazyList.continually(".").take(eightMB).mkString
    HashType.S3Etag.calculateHash(value) shouldBe "f89801f68b5028d64e0238ffb5a1b8e0"
  }

  it should "calculate an etag hash on long data" in {
    val eightMB = 8 * 1024 * 1024
    val value = LazyList.continually(".").take(eightMB + 1).mkString
    HashType.S3Etag.calculateHash(value) shouldBe "8e224b463f4f5202c9621820f7690a01-2"
  }

  it should "calculate an md5 hash" in {
    HashType.Md5.calculateHash("hello") shouldBe "5d41402abc4b2a76b9719d911017c592"
  }

  it should "calculate a sha256 hash" in {
    HashType.Sha256.calculateHash("hello") shouldBe "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824"
  }
}
