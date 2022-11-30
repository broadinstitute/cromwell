package cromwell.filesystems.drs

import cloud.nio.impl.drs.{MockDrsCloudNioFileSystemProvider, MockDrsPaths}
import com.typesafe.config.{Config, ConfigFactory}
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.jdk.CollectionConverters._


class DrsResolverSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  private val drsResolverConfig: Config = ConfigFactory.parseMap(
    Map(
      "resolver.url" -> "https://drshub-url/drshub_v4",
      "access-token-acceptable-ttl" -> "1 minute"
    ).asJava
  )

  private val mockFileSystemProvider = new MockDrsCloudNioFileSystemProvider(config = drsResolverConfig)
  private val drsPathBuilder = DrsPathBuilder(mockFileSystemProvider, None)


  behavior of "DrsResolver"

  it should "find DRS path from a GCS path" in {
    val drsPath = drsPathBuilder.build(MockDrsPaths.drsPathResolvingGcsPath).get.asInstanceOf[DrsPath]

    DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync() should be (MockDrsPaths.drsRelativePath)
  }

  it should "find DRS path from a path replacing characters" in {
    val drsPath = drsPathBuilder.build(MockDrsPaths.drsPathWithNonPathChars).get.asInstanceOf[DrsPath]

    DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync() should be (MockDrsPaths.drsReplacedChar)
  }

  it should "find DRS path from a file name" in {
    val drsPath = drsPathBuilder.build(MockDrsPaths.drsPathResolvingWithFileName).get.asInstanceOf[DrsPath]

    DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync() should be (MockDrsPaths.gcsRelativePathWithFileName)
  }

  it should "find DRS path from a localization path" in {
    val drsPath = drsPathBuilder.build(MockDrsPaths.drsPathResolvingWithLocalizationPath).get.asInstanceOf[DrsPath]

    DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync() should be (MockDrsPaths.gcsRelativePathWithFileNameFromLocalizationPath)
  }

  it should "find DRS path from all the paths" in {
    val drsPath = drsPathBuilder.build(MockDrsPaths.drsPathResolvingWithAllThePaths).get.asInstanceOf[DrsPath]

    DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync() should be (MockDrsPaths.gcsRelativePathWithFileNameFromAllThePaths)
  }

  it should "throw GcsUrlNotFoundException when DRS path doesn't resolve to at least one GCS url" in {
    val drsPath = drsPathBuilder.build(MockDrsPaths.drsPathResolvingToNoGcsPath).get.asInstanceOf[DrsPath]

    the[RuntimeException] thrownBy {
      DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync()
    } should have message s"Error while resolving DRS path: ${drsPath.pathAsString}. Error: UrlNotFoundException: No gs url associated with given DRS path."
  }

  it should "throw Runtime Exception when the DRS Resolver can't find the given DRS path " in {
    val drsPath = drsPathBuilder.build(MockDrsPaths.drsPathNotExistingInDrsResolver).get.asInstanceOf[DrsPath]

    the[RuntimeException] thrownBy {
      DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync()
    } should have message
      s"Error while resolving DRS path: ${drsPath.pathAsString}. " +
        s"Error: RuntimeException: Unexpected response resolving ${drsPath.pathAsString} " +
        s"through DRS Resolver url https://drshub-url/drshub_v4. Error: 404 Not Found."
  }
}
