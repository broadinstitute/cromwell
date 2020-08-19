package cloud.nio.impl.ftp

import java.io.IOException

import cats.effect.IO
import cloud.nio.impl.ftp.FtpUtil._
import org.apache.commons.net.ftp.FTPClient
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.specs2.mock.Mockito

import scala.concurrent.duration._


class FtpUtilSpec extends AnyFlatSpec with Matchers with Mockito {

  behavior of "autoRelease"

  it should "release the lease when the client fails the operation without throwing" in {
    val clientPool = new FtpClientPool(1, 10.minutes, () => { new FTPClient })
    val lease = clientPool.acquire()
    
    val action = autoRelease(IO.pure(lease)) { _ =>
      IO.raiseError(FtpIoException("boom", 1, "re-boom"))
    }
    an[FtpIoException] shouldBe thrownBy(action.unsafeRunSync())
    an[IllegalStateException] shouldBe thrownBy(lease.get())
    clientPool.live() shouldBe 1
  }

  it should "invalidate the lease when the client fails the operation by throwing" in {
    val clientPool = new FtpClientPool(1, 10.minutes, () => { new FTPClient })
    val lease = clientPool.acquire()
    
    val action = autoRelease(IO.pure(lease)) { _ =>
      IO.raiseError(FtpIoException("boom", 1, "re-boom", Option(new IOException("baaaam"))))
    }
    an[FtpIoException] shouldBe thrownBy(action.unsafeRunSync())
    an[IllegalStateException] shouldBe thrownBy(lease.get())
    clientPool.live() shouldBe 0
  }

  it should "release the lease when the operation succeeds" in {
    val clientPool = new FtpClientPool(1, 10.minutes, () => { new FTPClient })
    val lease = clientPool.acquire()
    
    val action = autoRelease(IO.pure(lease)) { _ =>
      IO.pure("yeahh")
    }
    action.unsafeRunSync() shouldBe "yeahh"
    an[IllegalStateException] shouldBe thrownBy(lease.get())
    clientPool.live() shouldBe 1
  }
}
