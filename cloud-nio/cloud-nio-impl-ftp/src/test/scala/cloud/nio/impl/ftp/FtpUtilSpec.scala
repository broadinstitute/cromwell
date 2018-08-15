package cloud.nio.impl.ftp

import java.nio.file.{FileAlreadyExistsException, NoSuchFileException}

import cloud.nio.impl.ftp.FtpUtil.{FtpIoException, FtpOperation}
import org.apache.commons.net.ftp.{FTPClient, FTPReply}
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito
import scala.concurrent.duration._

class FtpUtilSpec extends FlatSpec with Matchers with Mockito {

  behavior of "FtpUtilSpec"

  it should "generate somewhat accurate exceptions" in {
    val client = mock[FTPClient]
    val operation = FtpOperation("ftp.example.com", "location", "do something")
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
  
  it should "release the lease when the client fails the operation without throwing" in {
    val clientPool = new FtpClientPool(1, 10.minutes, () => { new FTPClient })
    val operation = FtpOperation("ftp.example.com", "location", "do something")
    val lease = clientPool.acquire()
    an[FtpIoException] shouldBe thrownBy(operation.fail(lease, None))
    an[IllegalStateException] shouldBe thrownBy(lease.get())
    clientPool.live() shouldBe 1
  }

  it should "invalidate the lease when the client fails the operation by throwing" in {
    val clientPool = new FtpClientPool(1, 10.minutes, () => { new FTPClient })
    val operation = FtpOperation("ftp.example.com", "location", "do something")
    val lease = clientPool.acquire()
    an[FtpIoException] shouldBe thrownBy(operation.fail(lease, Option(new Exception)))
    an[IllegalStateException] shouldBe thrownBy(lease.get())
    clientPool.live() shouldBe 0
  }

}
