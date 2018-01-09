package cromwell.filesystems.oss.nio

import cromwell.core.TestKitSuite

class OssStorageFileSystemSpec extends TestKitSuite with OssNioUtilSpec {
  behavior of s"OssStorageFileSystemSpec"

  it should "get right path" in {
    val ossPath = mockFileSystem.getPath("/test-file-system")
    ossPath.bucket shouldEqual(bucket)
    ossPath.key shouldEqual("test-file-system")
  }

  it should "has right view name" in {
    val fs = mockFileSystem

    fs.supportedFileAttributeViews should contain (OssStorageFileSystem.BASIC_VIEW)
    fs.supportedFileAttributeViews should contain (OssStorageFileSystem.OSS_VIEW)
  }

  it should "do not support some method" in {
    an [UnsupportedOperationException] should be thrownBy mockFileSystem.newWatchService
  }

  it should "return some expected simple mocked result" in {
    mockFileSystem.isOpen shouldBe true
    mockFileSystem.isReadOnly shouldBe false
    mockFileSystem.getSeparator shouldBe OssStorageFileSystem.SEPARATOR
  }
}
