package cloud.nio.impl.ftp.operations

import java.nio.file.{FileAlreadyExistsException, NoSuchFileException}

import cloud.nio.impl.ftp.FtpUtil.FtpIoException
import common.assertion.CromwellTimeoutSpec
import org.apache.commons.net.ftp.{FTPClient, FTPReply}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.specs2.mock.Mockito

class FtpOperationSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with Mockito {
  it should "generate somewhat accurate exceptions" in {
    val client = mock[FTPClient]
    val operation = FtpListFiles("ftp.example.com", "location", "do something")
    client.getReplyCode.returns(FTPReply.FILE_UNAVAILABLE)

    client.getReplyString.returns("no such file")
    operation.generateException(client, None) shouldBe a[NoSuchFileException]

    client.getReplyString.returns("file already exists")
    operation.generateException(client, None) shouldBe a[FileAlreadyExistsException]

    client.getReplyString.returns("blah")
    val withoutCause = operation.generateException(client, None)
    withoutCause shouldBe an[FtpIoException]
    withoutCause.getCause shouldBe null

    val cause = new Exception()
    val withCause = operation.generateException(client, Option(cause))
    withCause shouldBe an[FtpIoException]
    withCause.getCause shouldBe cause
  }
}
