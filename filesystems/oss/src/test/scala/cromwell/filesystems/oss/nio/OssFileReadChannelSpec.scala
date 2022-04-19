package cromwell.filesystems.oss.nio

import java.nio.charset.Charset

import cromwell.core.TestKitSuite
import org.scalatest.{BeforeAndAfter}

import scala.util.Try
import scala.util.control.Breaks

object OssFileReadChannelSpec {
  val FILENAME = "/test-oss-read-file"
  val CONTENT = "Hello World!"

  implicit class Crossable[X](xs: Iterable[X]) {
    def cross[Y](ys: Iterable[Y]) = for { x <- xs; y <- ys } yield (x, y)
  }
}

class OssFileReadChannelSpec extends TestKitSuite with OssNioUtilSpec with BeforeAndAfter {
  behavior of s"OssFileReadChannelSpec"

  import OssFileReadChannelSpec._


  def getPath = OssStoragePath.getPath(ossFileSystem, FILENAME)

  before {
     Try(OssAppendOutputStream(ossClient, getPath, true)) foreach {_.write(CONTENT.getBytes("UTF-8"))}
  }

  after {
    Try(deleteObject(getPath))
  }

  it should "has the right size" taggedAs NeedAK in {
    val channel = OssFileReadChannel(ossClient, 0L, getPath)
    channel.size shouldEqual(CONTENT.length)
  }

  it should "has the right content" taggedAs NeedAK in {
    List.range(1, CONTENT.length + 1) foreach { bufferSize =>verifySameContent(bufferSize)}
    for (bufferSize <- List.range(1, CONTENT.length + 1); position <- List.range(0, CONTENT.length)) {
      verifySameContent(bufferSize, position.toLong)
    }
  }

  it should "has the right position after seeking" taggedAs NeedAK in {
    val channel = OssFileReadChannel(ossClient, 0L, getPath)
    channel.size shouldEqual(CONTENT.length)

    channel.position(1)

    channel.position shouldEqual(1)
  }

  def verifySameContent(bufferSize: Int, position: Long = 0) = {
    val channel = OssFileReadChannel(ossClient, position, getPath)

    import java.nio.ByteBuffer
    val buf = ByteBuffer.allocate(bufferSize)

    val loop = new Breaks
    val builder = new StringBuilder

    var bytesRead = channel.read(buf)
    loop.breakable {
      while (bytesRead != -1) {
        buf.flip()
        val charset = Charset.forName("UTF-8");

        builder.append(charset.decode(buf).toString())
        buf.clear
        bytesRead = channel.read(buf)
      }
    }

    builder.toString shouldEqual CONTENT.substring(position.toInt)
  }
}
