package cromwell.filesystems.drs

import com.typesafe.config.{Config, ConfigFactory}
import org.scalatest.{FlatSpec, Matchers}


class DrsResolverSpec extends FlatSpec with Matchers {

  private val marthaConfig: Config = ConfigFactory.parseString(
    """martha {
      |   url = "http://matha-url"
      |   request.json-template = "{"key": "${holder}"}"
      |}
      |""".stripMargin
  )

  private val mockFileSystemProvider = new MockDrsCloudNioFileSystemProvider(marthaConfig)
  private val drsPathBuilder = DrsPathBuilder(mockFileSystemProvider)

  val gcsRelativePath = "mybucket/foo.txt"


  behavior of "DrsResolver"

  it should "find GCS path when its the only one in url array" in {
    val drsPath = drsPathBuilder.build(MockDrsPaths.drsPathResolvingToOneGcsPath).get.asInstanceOf[DrsPath]

    DrsResolver.getContainerRelativePath(drsPath) should be (gcsRelativePath)
  }

  it should "find GCS path when DRS path resolves to multiple urls" in {
    val drsPath = drsPathBuilder.build(MockDrsPaths.drsPathResolvingToMultiplePaths).get.asInstanceOf[DrsPath]

    DrsResolver.getContainerRelativePath(drsPath) should be (gcsRelativePath)
  }

  it should "throw GcsUrlNotFoundException when DRS path doesn't resolve to at least one GCS url" in {
    val drsPath = drsPathBuilder.build(MockDrsPaths.drsPathResolvingToNoGcsPath).get.asInstanceOf[DrsPath]

    the[GcsUrlNotFoundException] thrownBy {
      DrsResolver.getContainerRelativePath(drsPath)
    } should have message s"DRS was not able to find a gs url associated with ${drsPath.pathAsString}."
  }

  it should "throw Runtime Exception when Martha can't find the given DRS path " in {
    val drsPath = drsPathBuilder.build(MockDrsPaths.drsPathNotExistingInMartha).get.asInstanceOf[DrsPath]

    the[RuntimeException] thrownBy {
      DrsResolver.getContainerRelativePath(drsPath)
    } should have message s"Unexpected response resolving ${drsPath.pathAsString} through Martha url http://matha-url. Error: 502 Bad Gateway."
  }
}
