package cromwell.filesystems.oss

import cromwell.core.TestKitSuite
import cromwell.filesystems.oss.nio.OssNioUtilSpec
import org.scalatest.{BeforeAndAfter, FlatSpecLike, Matchers}
import org.scalatest.TryValues._

class OssPathBuilderSpec extends TestKitSuite with FlatSpecLike with Matchers with OssNioUtilSpec with BeforeAndAfter {

  behavior of s"OssPathBuilerSpec"
  val testPathBuiler = OssPathBuilder(mockOssConf)

  it should "throw when no bucket in URI" in {
    testPathBuiler.build("oss:").failed.get shouldBe an[IllegalArgumentException]
    testPathBuiler.build("oss://").failed.get shouldBe an[IllegalArgumentException]
  }

  it should "throw when path has an invalid schema" in {
    testPathBuiler.build(s"gcs://$bucket$fileName").failed.get shouldBe an[IllegalArgumentException]
  }

  it should "has an empty key when no path specified" in {
    testPathBuiler.build(s"oss://$bucket").success.value.bucket shouldBe bucket
    testPathBuiler.build(s"oss://$bucket").success.value.key shouldBe empty
  }

  it should "start with separator when path specified" in {
    val path = testPathBuiler.build(s"oss://$bucket$fileName").success.value
    path.bucket shouldBe bucket
    path.nioPath.toString shouldBe fileName
    path.key shouldBe fileName.stripPrefix("/")
    path.pathAsString shouldBe s"oss://$bucket$fileName"
    path.pathWithoutScheme shouldBe s"$bucket$fileName"
  }
 }
