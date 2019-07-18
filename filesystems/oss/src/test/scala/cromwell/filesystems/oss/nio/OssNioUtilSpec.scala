package cromwell.filesystems.oss.nio

import java.io.ByteArrayInputStream

import com.aliyun.oss.OSSClient
import org.scalatest._
import org.scalatest.mockito.MockitoSugar

import scala.util.control.Breaks
import scala.util.{Failure, Success, Try}

object NeedAK extends Tag("this test need oss storage access id and key")

object OssNioUtilSpec {
  val DEFAULT_BUCKET = "bcs-bucket"

  val DEFAULT_FILE_NAME = "/bcs-dir/bcs-file"

  val DEFAULT_CONTENT = "Hello World!"

  val ossInfo: Map[String, String] = Map(
    "endpoint" -> "",
    "access-id" -> "",
    "access-key" -> "",
    "bucket" -> DEFAULT_BUCKET
  )
}

trait OssNioUtilSpec extends FlatSpecLike with MockitoSugar with Matchers {

  override def withFixture(test: NoArgTest): Outcome = {
    if (test.tags.contains(NeedAK.name)) {
      Try(ossConf) match {
        case Success(_) => super.withFixture(test)
        case Failure(_) => cancel(NeedAK.name)
      }
    } else {
      super.withFixture(test)
    }
  }

  import OssNioUtilSpec._

  lazy val bucket: String = {
    val bucket = ossInfo.getOrElse("bucket", "mock-bucket")
    if (bucket.isEmpty) {
      throw new IllegalArgumentException("test bucket can not be empty")
    }

    bucket
  }

  lazy val ossConf: OssStorageConfiguration = Try{
    OssStorageConfiguration.parseMap(ossInfo)
  } getOrElse(throw new IllegalArgumentException("you should supply oss info before testing oss related operation"))

  lazy val mockOssConf: OssStorageConfiguration = new DefaultOssStorageConfiguration("mock-endpoint", "mock-id", "mock-key", None)

  lazy val ossProvider = OssStorageFileSystemProvider(ossConf)
  lazy val mockProvider = OssStorageFileSystemProvider(mockOssConf)
  lazy val ossFileSystem = OssStorageFileSystem(bucket, ossConf)
  lazy val mockFileSystem = OssStorageFileSystem(bucket, mockOssConf)
  val fileName = DEFAULT_FILE_NAME
  val fileContent = DEFAULT_CONTENT

  lazy val ossClient: OSSClient = mockOssConf.newOssClient()

  def contentAsString(path: OssStoragePath): String = {
    val ossObject = ossClient.getObject(path.bucket, path.key)

    val in = OssStorageRetry.from(
      () => ossObject.getObjectContent
    )

    val maxLen = 1024
    val loop = new Breaks
    val result = new StringBuilder
    loop.breakable {
      while(true) {
        val b = new Array[Byte](maxLen)
        val got = OssStorageRetry.from(
          () => in.read(b, 0, maxLen)
        )
        if (got <= 0) {
          loop.break()
        }
        result.append(new String(b, 0, got))
      }
    }

    result.toString()
  }

  def deleteObject(path: OssStoragePath): Unit = {
    OssStorageRetry.from(
      () => ossClient.deleteObject(path.bucket, path.key)
    )
  }

  def writeObject(path: OssStoragePath): Unit = {
    OssStorageRetry.from{
      () => ossClient.putObject(path.bucket, path.key, new ByteArrayInputStream(fileContent.getBytes()))
    }
    ()
  }
}
