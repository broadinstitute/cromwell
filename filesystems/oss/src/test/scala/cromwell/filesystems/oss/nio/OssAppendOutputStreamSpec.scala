package cromwell.filesystems.oss.nio

import cromwell.core.TestKitSuite

class OssAppendOutputStreamSpec extends TestKitSuite with OssNioUtilSpec  {

  behavior of s"OssAppendOutputStream"

  "write batch" should "work" taggedAs NeedAK in {
    val path = OssStoragePath.getPath(ossFileSystem, "/test-oss-append")
    val stream = OssAppendOutputStream(ossClient, path, true)

    val content: String = "haha"
    stream.write(content.getBytes)

    contentAsString(path) shouldEqual content
    stream.position shouldEqual content.length
  }

  "write single" should "work" taggedAs NeedAK in {
    val c: Char = 'c'
    val path = OssStoragePath.getPath(ossFileSystem, "/test-oss-append")
    val stream = OssAppendOutputStream(ossClient, path, true)

    stream.write(c.toInt)

    contentAsString(path) shouldEqual c.toString
    stream.position shouldEqual 1
  }

  "write range" should "work" taggedAs NeedAK in {
    val path = OssStoragePath.getPath(ossFileSystem, "/test-oss-append")
    val stream = OssAppendOutputStream(ossClient, path, true)

    val content: String = "haha"
    stream.write(content.getBytes, 1, 1)

    contentAsString(path) shouldEqual 'a'.toString
    stream.position shouldEqual 1
  }

}
