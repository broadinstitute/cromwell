package cromwell.backend.impl.sfs.config

import java.util.UUID

import akka.event.LoggingAdapter
import better.files._
import com.typesafe.config.{ConfigFactory, ConfigValueFactory}
import cromwell.backend.callcaching.FileHashingActor.SingleFileHashRequest
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.sfs.SharedFileSystemBackendInitializationData
import cromwell.core.path.DefaultPathBuilder
import org.apache.commons.codec.digest.DigestUtils
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import org.specs2.mock.Mockito
import wdl4s.values.WdlFile

import scala.util.Success

class ConfigHashingStrategySpec extends FlatSpec with Matchers with TableDrivenPropertyChecks with Mockito with BeforeAndAfterAll {

  behavior of "ConfigHashingStrategy"

  val steak = "Steak"
  val steakHash = DigestUtils.md5Hex(steak)
  val file = File.newTemporaryFile()
  val symLinksDir = File.newTemporaryDirectory("sym-dir")
  val pathHash = DigestUtils.md5Hex(file.pathAsString)
  val md5File = file.sibling(s"${file.name}.md5")
  // Not the md5 value of "Steak". This is intentional so we can verify which hash is used depending on the strategy
  val md5FileHash = "103508832bace55730c8ee8d89c1a45f"

  override def beforeAll() = {
    file.write(steak)
    ()
  }

  private def randomName(): String = UUID.randomUUID().toString

  def mockRequest(withSibling: Boolean, symlink: Boolean = false) = {
    if (withSibling && md5File.notExists) md5File.write(md5FileHash)
    val request = mock[SingleFileHashRequest]
    val requestFile = if (symlink) {
      val symLink : File = symLinksDir./(s"symlink-${randomName()}")
      symLink.symbolicLinkTo(file)
      symLink
    } else file

    val workflowPaths = mock[WorkflowPaths]
    workflowPaths.pathBuilders returns List(DefaultPathBuilder)

    val initData = mock[SharedFileSystemBackendInitializationData]
    initData.workflowPaths returns workflowPaths

    request.file returns WdlFile(requestFile.pathAsString)
    request.initializationData returns Option(initData)

    request
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
      (false, true, pathHash),
      (true, false, pathHash),
      (false, false, pathHash)
    )

    forAll(table) { (check, withMd5, expected) =>
      md5File.delete(swallowIOExceptions = true)
      val checkSibling = makeStrategy("path", Option(check))

      checkSibling.getHash(mockRequest(withMd5, symlink = false), mock[LoggingAdapter]) shouldBe Success(expected)

      val symLinkRequest: SingleFileHashRequest = mockRequest(withMd5, symlink = true)
      val symlink = File(symLinkRequest.file.valueString)

      symlink.isSymbolicLink shouldBe true
      DigestUtils.md5Hex(symlink.pathAsString) should not be expected
      checkSibling.getHash(symLinkRequest, mock[LoggingAdapter]) shouldBe Success(expected)
    }
  }

  it should "create a file hashing strategy from config" in {
    val defaultSibling = makeStrategy("file")
    defaultSibling.isInstanceOf[HashFileStrategy] shouldBe true
    defaultSibling.checkSiblingMd5 shouldBe false

    val checkSibling = makeStrategy("file", Option(true))

    checkSibling.isInstanceOf[HashFileStrategy] shouldBe true
    checkSibling.checkSiblingMd5 shouldBe true
    checkSibling.toString shouldBe "Call caching hashing strategy: Check first for sibling md5 and if not found hash file content."

    val dontCheckSibling = makeStrategy("file", Option(false))

    dontCheckSibling.isInstanceOf[HashFileStrategy] shouldBe true
    dontCheckSibling.checkSiblingMd5 shouldBe false
    dontCheckSibling.toString shouldBe "Call caching hashing strategy: hash file content."
  }

  it should "have a file hashing strategy and use md5 sibling file when appropriate" in {
    val table = Table(
      ("check", "withMd5", "expected"),
      (true, true, md5FileHash),
      (false, true, steakHash),
      (true, false, steakHash),
      (false, false, steakHash)
    )

    forAll(table) { (check, withMd5, expected) =>
      md5File.delete(swallowIOExceptions = true)
      val checkSibling = makeStrategy("file", Option(check))

      checkSibling.getHash(mockRequest(withMd5, symlink = false), mock[LoggingAdapter]) shouldBe Success(expected)

      val symLinkRequest: SingleFileHashRequest = mockRequest(withMd5, symlink = true)
      val symlink = File(symLinkRequest.file.valueString)

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
