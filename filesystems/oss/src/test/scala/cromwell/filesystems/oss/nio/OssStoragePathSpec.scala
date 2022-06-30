package cromwell.filesystems.oss.nio

import cromwell.core.TestKitSuite
import scala.jdk.CollectionConverters._


class OssStoragePathSpec extends TestKitSuite with OssNioUtilSpec {
  behavior of s"OssStoragePath"


  it should s"has the same bucket with file system" in {

    val path = OssStoragePath.getPath(mockFileSystem, fileName)

    path.bucket shouldBe bucket

    path.toAbsolutePath.toString shouldBe fileName
  }

  it should s"has a separator-removed key" in {
    val path = OssStoragePath.getPath(mockFileSystem, fileName)

    path.key shouldBe fileName.stripPrefix(UnixPath.SEPARATOR.toString)
  }

  "a not absolute oss path" should s"has a NullOssStoragePath root path" in {
    val path = OssStoragePath.getPath(mockFileSystem, fileName.stripPrefix(UnixPath.SEPARATOR.toString))

    path.getRoot shouldBe a [NullOssStoragePath]
  }

  "an absolute oss path" should s"has a OssStoragePathImpl root path" in {
    val path = OssStoragePath.getPath(mockFileSystem, fileName)

    path.getRoot shouldBe an [OssStoragePathImpl]
  }

  it should s"has right iterator" in {
    val path = OssStoragePath.getPath(mockFileSystem, fileName)

    var subs = List.empty[String]
    path.iterator().asScala foreach(p => subs = subs :+ p.toString)

    subs.head shouldBe "bcs-dir"
    subs(1) shouldBe "bcs-file"
  }

  it should s"has right relativize" in {
    val path = OssStoragePath.getPath(mockFileSystem, fileName)

    val path1 = OssStoragePath.getPath(mockFileSystem, "/bcs-dir/bcs-file1")

    path.relativize(path1).toString shouldEqual "../bcs-file1"

    val path2 = OssStoragePath.getPath(mockFileSystem, "/bcs-dir1/bcs-file2")
    path.relativize(path2).toString shouldEqual "../../bcs-dir1/bcs-file2"
  }

  it should s"has right pathAsString" in {
    val path = OssStoragePath.getPath(mockFileSystem, fileName)

    path.pathAsString shouldEqual s"oss://$bucket$fileName"
  }

}
