package cloud.nio.impl.ftp

import java.io.OutputStream

import common.assertion.CromwellTimeoutSpec
import org.apache.commons.net.ftp.FTPClient
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.specs2.mock.Mockito
import scala.concurrent.duration._


class LeaseOutputStreamSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with Mockito {

  behavior of "LeaseInputStreamSpec"

  it should "complete the command and release the lease when closing the stream" in {
    val os = new TestOutputStream
    val mockClient = mock[FTPClient]
    var completed: Boolean = false
    mockClient.completePendingCommand().returns({
      completed = true
      true
    })
    val clientPool = new FtpClientPool(1, 10.minutes, () => { mockClient })
    val lease = clientPool.acquire()
    val leasedOutputStream = new LeasedOutputStream("host", "path", os, lease)

    leasedOutputStream.close()

    completed shouldBe true
    os.isClosed shouldBe true
    // When accessing a released lease, get throws an IllegalStateException
    an[IllegalStateException] shouldBe thrownBy(lease.get())
  }
  
  private class TestOutputStream extends OutputStream {
    var closed = false
    override def close() = closed = true
    def isClosed = closed
    override def write(b: Int) = {}
  }
}
