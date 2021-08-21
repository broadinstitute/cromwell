package cromwell.filesystems.oss.nio

import java.nio.file.attribute.FileTime
import com.aliyun.oss.model.ObjectMetadata
import common.mock.MockSugar
import cromwell.core.TestKitSuite
import org.mockito.Mockito._

import java.text.SimpleDateFormat
import java.util.{Date, Locale}

object OssStorageObjectAttributesSpec extends MockSugar {
  val DEFAULT_BUCKET = "bcs-bucket"

  val DEFAULT_FILE_NAME = "/bcs-dir/bcs-file"

  val DEFAULT_LENGTH: Long = 2102784

  val DEFAULT_MODIFIED: Date = {
    val target = "Thu Dec 21 15:19:27 CST 2017"
    val df = new SimpleDateFormat("EEE MMM dd kk:mm:ss z yyyy", Locale.ENGLISH)
    df.parse(target)
  }

  val DEFAULT_ETAG = "F80066F040BDA4F991DB5F8AEC9905FB"

  val DEFAULT_CONTENT_DISPOSITION: String = null

  val DEFAULT_CACHE_CONTROL: String = null

  val DEFAULT_CONTENT_ENCODING: String = null

  val DEFAULT_CONTENT_TYPE = "application/x-msdownload"

  def getObjectMeta: ObjectMetadata = {
    val meta = mock[ObjectMetadata]

    when(meta.getContentDisposition).thenReturn(DEFAULT_CONTENT_DISPOSITION)
    when(meta.getContentEncoding).thenReturn(DEFAULT_CONTENT_ENCODING)
    when(meta.getCacheControl).thenReturn(DEFAULT_CACHE_CONTROL)
    when(meta.getLastModified).thenReturn(DEFAULT_MODIFIED)
    when(meta.getETag).thenReturn(DEFAULT_ETAG)
    when(meta.getContentType).thenReturn(DEFAULT_CONTENT_TYPE)
    when(meta.getContentLength).thenReturn(DEFAULT_LENGTH)
    when(meta.getExpirationTime).thenThrow(new NullPointerException())

    meta
  }
}

class OssStorageObjectAttributesSpec extends TestKitSuite with OssNioUtilSpec {

  behavior of s"OssStorageObjectAttributes"

  import OssStorageObjectAttributesSpec._

  def getObject: OssStoragePathImpl = {
    OssStoragePath.getPath(mockFileSystem, fileName)
  }

  def getDir: OssStoragePathImpl = {
    OssStoragePath.getPath(mockFileSystem, "/bcs-dir/")
  }

  "an oss object attr" should "be an right" in {
    val attr = OssStorageObjectAttributes(getObjectMeta, getObject)

    attr.fileKey shouldEqual getObject.pathAsString

    attr.creationTime shouldEqual attr.lastModifiedTime()
    attr.lastAccessTime shouldEqual FileTime.fromMillis(0)
    attr.cacheControl shouldBe empty
    attr.contentDisposition shouldBe empty
    attr.contentEncoding shouldBe empty
    attr.etag shouldBe Option(DEFAULT_ETAG)
    attr.size shouldBe DEFAULT_LENGTH
  }

  "an oss directory attr" should "be an right" in {
    val attr = OssStorageDirectoryAttributes(getDir)

    attr.fileKey shouldEqual getDir.pathAsString

    attr.creationTime shouldEqual attr.lastModifiedTime()
    attr.lastAccessTime shouldEqual FileTime.fromMillis(0)
    attr.cacheControl shouldBe empty
    attr.contentDisposition shouldBe empty
    attr.contentEncoding shouldBe empty
    attr.etag shouldBe empty
    attr.size shouldBe 0
  }
}
