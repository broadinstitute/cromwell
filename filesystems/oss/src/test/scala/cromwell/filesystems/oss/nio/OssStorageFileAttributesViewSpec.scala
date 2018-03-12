package cromwell.filesystems.oss.nio

import com.aliyun.oss.OSSClient
import org.mockito.Mockito._
import org.mockito.Matchers._


class OssStorageFileAttributesViewSpec extends OssNioUtilSpec {
  behavior of "OssStorageFileAttributesView"

  import OssStorageObjectAttributesSpec._

  def getObject = {
    OssStoragePath.getPath(mockFileSystem, fileName)
  }

  def getDir = {
    OssStoragePath.getPath(mockFileSystem, "/bcs-dir/")
  }

  it should "return an object attr" in {
    val ossClient = mock[OSSClient]
    when(ossClient.doesObjectExist(anyString(), anyString())).thenReturn(true)
    val meta = getObjectMeta
    when(ossClient.getObjectMetadata(anyString(), anyString())).thenReturn(meta)

    val view = OssStorageFileAttributesView(ossClient, getObject)
    view.readAttributes shouldBe an [OssStorageObjectAttributes]
  }

  it should "return an dir attr" in {
    val ossClient = mock[OSSClient]
    when(ossClient.doesObjectExist(anyString(), anyString())).thenReturn(true)
    val view = OssStorageFileAttributesView(ossClient, getDir)
    view.readAttributes shouldBe a [OssStorageDirectoryAttributes]
  }

}
