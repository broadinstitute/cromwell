package cromwell.backend.impl.sfs.config

import java.util.UUID

import akka.event.LoggingAdapter
import com.typesafe.config.{ConfigFactory, ConfigValueFactory}
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.standard.StandardInitializationData
import cromwell.backend.standard.callcaching.StandardFileHashingActor.SingleFileHashRequest
import cromwell.core.path.{DefaultPathBuilder, Path}
import org.apache.commons.codec.digest.DigestUtils
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import org.specs2.mock.Mockito
import wom.values.WomSingleFile
import scala.util.Success

class ConfigHashingStrategySpec extends FlatSpec with Matchers with TableDrivenPropertyChecks with Mockito with BeforeAndAfterAll {

  behavior of "ConfigHashingStrategy"

  val steak = "Steak"
  val steakMd5 = DigestUtils.md5Hex(steak)
  val steakXxh64 = HashFileXxH64StrategyMethods.xxh64sumString(steak)
  val file = DefaultPathBuilder.createTempFile()
  val symLinksDir = DefaultPathBuilder.createTempDirectory("sym-dir")
  val pathMd5 = DigestUtils.md5Hex(file.pathAsString)
  val md5File = file.sibling(s"${file.name}.md5")
  // Not the md5 value of "Steak". This is intentional so we can verify which hash is used depending on the strategy
  val md5FileHash = "103508832bace55730c8ee8d89c1a45f"

  override def beforeAll() = {
    file.write(steak)
    ()
  }

  private def randomName(): String = UUID.randomUUID().toString

  def mockRequest(withSibling: Boolean, symlink: Boolean) = {
    if (withSibling && md5File.notExists) md5File.write(md5FileHash + System.lineSeparator())
    val requestFile = if (symlink) {
      val symLink: Path = symLinksDir./(s"symlink-${randomName()}")
      symLink.symbolicLinkTo(file)
      symLink
    } else file

    val workflowPaths = mock[WorkflowPaths]
    workflowPaths.pathBuilders returns List(DefaultPathBuilder)

    val initData = mock[StandardInitializationData]
    initData.workflowPaths returns workflowPaths
    SingleFileHashRequest(null, null, WomSingleFile(requestFile.pathAsString), Option(initData))
  }

  def makeStrategy(strategy: String, checkSibling: Option[Boolean] = None) = {
    val conf = ConfigFactory.parseString(s"""hashing-strategy: "$strategy"""")
    ConfigHashingStrategy(
      checkSibling map { check => conf.withValue("check-sibling-md5", ConfigValueFactory.fromAnyRef(check)) } getOrElse conf
    )
  }

  it should "create a path hashing strategy from config" in {
    val defaultSibling = makeStrategy("path")
    defaultSibling.isInstanceOf[HashPathStrategy] shouldBe true
    defaultSibling.checkSiblingMd5 shouldBe false

    val checkSibling = makeStrategy("path", Option(true))

    checkSibling.isInstanceOf[HashPathStrategy] shouldBe true
    checkSibling.checkSiblingMd5 shouldBe true
    checkSibling.toString shouldBe "Call caching hashing strategy: Check first for sibling md5 and if not found hash file path."

    val dontCheckSibling = makeStrategy("path", Option(false))

    dontCheckSibling.isInstanceOf[HashPathStrategy] shouldBe true
    dontCheckSibling.checkSiblingMd5 shouldBe false
    dontCheckSibling.toString shouldBe "Call caching hashing strategy: hash file path."
  }

  it should "have a path hashing strategy and use md5 sibling file when appropriate" in {
    val table = Table(
      ("check", "withMd5", "expected"),
      (true, true, md5FileHash),
      (false, true, pathMd5),
      (true, false, pathMd5),
      (false, false, pathMd5)
    )

    forAll(table) { (check, withMd5, expected) =>
      md5File.delete(swallowIOExceptions = true)
      val checkSibling = makeStrategy("path", Option(check))

      checkSibling.getHash(mockRequest(withMd5, symlink = false), mock[LoggingAdapter]) shouldBe Success(expected)

      val symLinkRequest: SingleFileHashRequest = mockRequest(withMd5, symlink = true)
      val symlink = DefaultPathBuilder.get(symLinkRequest.file.valueString)

      symlink.isSymbolicLink shouldBe true
      DigestUtils.md5Hex(symlink.pathAsString) should not be expected
      checkSibling.getHash(symLinkRequest, mock[LoggingAdapter]) shouldBe Success(expected)
    }
  }

  it should "create a path+modtime hashing strategy from config" in {
    val defaultSibling = makeStrategy("path+modtime")
    defaultSibling.isInstanceOf[HashPathModTimeStrategy] shouldBe true
    defaultSibling.checkSiblingMd5 shouldBe false

    val checkSibling = makeStrategy("path+modtime", Option(true))

    checkSibling.isInstanceOf[HashPathModTimeStrategy] shouldBe true
    checkSibling.checkSiblingMd5 shouldBe true
    checkSibling.toString shouldBe "Call caching hashing strategy: Check first for sibling md5 and if not found hash file path and last modified time."

    val dontCheckSibling = makeStrategy("path+modtime", Option(false))

    dontCheckSibling.isInstanceOf[HashPathModTimeStrategy] shouldBe true
    dontCheckSibling.checkSiblingMd5 shouldBe false
    dontCheckSibling.toString shouldBe "Call caching hashing strategy: hash file path and last modified time."
  }

  it should "have a path+modtime hashing strategy and use md5 sibling file when appropriate" in {
    // Have to define this here to make sure the timestamp is correct. Since the beforeAll() function modifies the file.
    val pathModTimeHash = DigestUtils.md5Hex(file.pathAsString + file.lastModifiedTime.toString)
    val table = Table(
      ("check", "withMd5", "expected"),
      (true, true, md5FileHash),
      (false, true, pathModTimeHash),
      (true, false, pathModTimeHash),
      (false, false, pathModTimeHash)
    )

    forAll(table) { (check, withMd5, expected) =>
      md5File.delete(swallowIOExceptions = true)
      val checkSibling = makeStrategy("path+modtime", Option(check))

      checkSibling.getHash(mockRequest(withMd5, symlink = false), mock[LoggingAdapter]) shouldBe Success(expected)

      val symLinkRequest: SingleFileHashRequest = mockRequest(withMd5, symlink = true)
      val symlink = DefaultPathBuilder.get(symLinkRequest.file.valueString)

      symlink.isSymbolicLink shouldBe true
      DigestUtils.md5Hex(symlink.pathAsString) should not be expected
      checkSibling.getHash(symLinkRequest, mock[LoggingAdapter]) shouldBe Success(expected)
    }
  }

  it should "create a md5 hashing strategy from config" in {
    val defaultSibling = makeStrategy("file")
    defaultSibling.isInstanceOf[HashFileMd5Strategy] shouldBe true
    defaultSibling.checkSiblingMd5 shouldBe false

    val checkSibling = makeStrategy("md5", Option(true))

    checkSibling.isInstanceOf[HashFileMd5Strategy] shouldBe true
    checkSibling.checkSiblingMd5 shouldBe true
    checkSibling.toString shouldBe "Call caching hashing strategy: Check first for sibling md5 and if not found hash file content with md5."

    val dontCheckSibling = makeStrategy("file", Option(false))

    dontCheckSibling.isInstanceOf[HashFileMd5Strategy] shouldBe true
    dontCheckSibling.checkSiblingMd5 shouldBe false
    dontCheckSibling.toString shouldBe "Call caching hashing strategy: hash file content with md5."
  }

  it should "have a file hashing strategy and use md5 sibling file when appropriate" in {
    val table = Table(
      ("check", "withMd5", "expected"),
      (true, true, md5FileHash),
      (false, true, steakMd5),
      (true, false, steakMd5),
      (false, false, steakMd5)
    )

    forAll(table) { (check, withMd5, expected) =>
      md5File.delete(swallowIOExceptions = true)
      val checkSibling = makeStrategy("file", Option(check))

      checkSibling.getHash(mockRequest(withMd5, symlink = false), mock[LoggingAdapter]) shouldBe Success(expected)

      val symLinkRequest: SingleFileHashRequest = mockRequest(withMd5, symlink = true)
      val symlink = DefaultPathBuilder.get(symLinkRequest.file.valueString)

      symlink.isSymbolicLink shouldBe true
      checkSibling.getHash(symLinkRequest, mock[LoggingAdapter]) shouldBe Success(expected)
    }
  }

  it should "create a xxh64 hashing strategy from config" in {
    val defaultSibling = makeStrategy("xxh64")
    defaultSibling.isInstanceOf[HashFileXxH64Strategy] shouldBe true
    defaultSibling.checkSiblingMd5 shouldBe false

    val checkSibling = makeStrategy("xxh64", Option(true))

    checkSibling.isInstanceOf[HashFileXxH64Strategy] shouldBe true
    checkSibling.checkSiblingMd5 shouldBe true
    checkSibling.toString shouldBe "Call caching hashing strategy: Check first for sibling md5 and if not found hash file content with xxh64."

    val dontCheckSibling = makeStrategy("xxh64", Option(false))

    dontCheckSibling.isInstanceOf[HashFileXxH64Strategy] shouldBe true
    dontCheckSibling.checkSiblingMd5 shouldBe false
    dontCheckSibling.toString shouldBe "Call caching hashing strategy: hash file content with xxh64."
  }

  it should "have a xxh64 hashing strategy and use md5 sibling file when appropriate" in {
    val table = Table(
      ("check", "withMd5", "expected"),
      (true, true, md5FileHash),
      (false, true, steakXxh64),
      (true, false, steakXxh64),
      (false, false, steakXxh64)
    )

    forAll(table) { (check, withMd5, expected) =>
      md5File.delete(swallowIOExceptions = true)
      val checkSibling = makeStrategy("xxh64", Option(check))

      checkSibling.getHash(mockRequest(withMd5, symlink = false), mock[LoggingAdapter]) shouldBe Success(expected)

      val symLinkRequest: SingleFileHashRequest = mockRequest(withMd5, symlink = true)
      val symlink = DefaultPathBuilder.get(symLinkRequest.file.valueString)

      symlink.isSymbolicLink shouldBe true
      checkSibling.getHash(symLinkRequest, mock[LoggingAdapter]) shouldBe Success(expected)
    }
  }

  it should "create a fingerprint strategy from config" in {
    val defaultFingerprint:  FingerprintStrategy = makeStrategy("fingerprint").asInstanceOf[FingerprintStrategy]
    defaultFingerprint.isInstanceOf[FingerprintStrategy] shouldBe true
    defaultFingerprint.checkSiblingMd5 shouldBe false
    defaultFingerprint.fingerprintSize shouldBe 10 * 1024 * 1024

    val config = ConfigFactory.parseString(
      """|hashing-strategy: "fingerprint"
         |fingerprint-size: 123456789
         |""".stripMargin)
    val otherFingerprint: FingerprintStrategy = ConfigHashingStrategy.apply(config).asInstanceOf[FingerprintStrategy]
    otherFingerprint.fingerprintSize shouldBe 123456789
    otherFingerprint.isInstanceOf[FingerprintStrategy] shouldBe true

    val checkSibling = makeStrategy("fingerprint", Option(true))

    checkSibling.isInstanceOf[FingerprintStrategy] shouldBe true
    checkSibling.checkSiblingMd5 shouldBe true
    checkSibling.toString shouldBe "Call caching hashing strategy: Check first for sibling md5 and if not found fingerprint the file with last modified time, size and a xxh64 hash of the first part of the file."

    val dontCheckSibling = makeStrategy("fingerprint", Option(false))

    dontCheckSibling.isInstanceOf[FingerprintStrategy] shouldBe true
    dontCheckSibling.checkSiblingMd5 shouldBe false
    dontCheckSibling.toString shouldBe "Call caching hashing strategy: fingerprint the file with last modified time, size and a xxh64 hash of the first part of the file."

  }

  it should "have a fingerprint strategy and use md5 sibling file when appropriate" in {
    val fingerPrintHash = HashFileXxH64StrategyMethods.xxh64sumString(file.lastModifiedTime.toEpochMilli.toHexString +
                                                              file.size.toHexString) + steakXxh64
    val table = Table(
      ("check", "withMd5", "expected"),
      (true, true, md5FileHash),
      (false, true, fingerPrintHash),
      (true, false, fingerPrintHash),
      (false, false, fingerPrintHash)
    )

    forAll(table) { (check, withMd5, expected) =>
      md5File.delete(swallowIOExceptions = true)
      val checkSibling = makeStrategy("fingerprint", Option(check))

      checkSibling.getHash(mockRequest(withMd5, symlink = false), mock[LoggingAdapter]) shouldBe Success(expected)

      val symLinkRequest: SingleFileHashRequest = mockRequest(withMd5, symlink = true)
      val symlink = DefaultPathBuilder.get(symLinkRequest.file.valueString)

      symlink.isSymbolicLink shouldBe true
      checkSibling.getHash(symLinkRequest, mock[LoggingAdapter]) shouldBe Success(expected)
    }
  }

  override def afterAll() = {
    file.delete(true)
    md5File.delete(true)
    ()
  }
}
