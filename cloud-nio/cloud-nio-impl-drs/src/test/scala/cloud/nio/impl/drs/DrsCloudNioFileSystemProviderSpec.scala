package cloud.nio.impl.drs

import org.scalatest.matchers.should.Matchers

class DrsCloudNioFileSystemProviderSpec extends org.scalatest.flatspec.AnyFlatSpec with Matchers {

  behavior of "DrsCloudNioFileSystemProvider"

  it should "checking is a directory exists should return false" in {
    val fileSystemProvider = new MockDrsCloudNioFileSystemProvider()
    val path = fileSystemProvider.getCloudNioPath("drs://foo/bar/")
    fileSystemProvider.checkDirectoryExists(path) should be(false)
  }

  it should "not delete" in {
    val fileSystemProvider = new MockDrsCloudNioFileSystemProvider()
    val path = fileSystemProvider.getCloudNioPath("drs://foo/bar/")
    the[UnsupportedOperationException] thrownBy {
      fileSystemProvider.deleteIfExists(path)
    } should have message("DRS currently doesn't support delete.")
  }
}
