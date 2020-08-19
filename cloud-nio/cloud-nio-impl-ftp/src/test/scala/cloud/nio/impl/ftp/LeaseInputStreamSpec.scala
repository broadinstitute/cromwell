package cloud.nio.impl.ftp

import java.io.InputStream

import org.apache.commons.net.ftp.FTPClient
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.specs2.mock.Mockito

import scala.concurrent.duration._


class LeaseInputStreamSpec extends AnyFlatSpec with Matchers with Mockito {

  behavior of "LeaseInputStreamSpec"

  it should "complete the command and release the lease when closing the stream" in {
    val is = new InputStream {
      var counter = 1
      override def read() = 1
      override def close() = counter = 0
      override def available() = counter
    }
    val mockClient = mock[FTPClient]
    var completed: Boolean = false
    mockClient.completePendingCommand().returns({
      completed = true
      true
    })
    val clientPool = new FtpClientPool(1, 10.minutes, () => { mockClient })
    val lease = clientPool.acquire()
    val leasedInputStream = new LeasedInputStream("host", "path", is, lease)

    leasedInputStream.close()

    completed shouldBe true
    is.available() shouldBe 0
    // When accessing a released lease, get throws an IllegalStateException
    an[IllegalStateException] shouldBe thrownBy(lease.get())
  }
}
