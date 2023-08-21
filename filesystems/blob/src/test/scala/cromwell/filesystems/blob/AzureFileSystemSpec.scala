package cromwell.filesystems.blob

import com.azure.storage.blob.nio.{AzureFileSystem, AzureFileSystemProvider}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.time.Instant
import scala.compat.java8.OptionConverters._
import scala.jdk.CollectionConverters._

class AzureFileSystemSpec extends AnyFlatSpec with Matchers {
  val now = Instant.now()
  val container = BlobContainerName("testConainer")
  val exampleSas = BlobPathBuilderFactorySpec.buildExampleSasToken(now)
  val exampleConfig = BlobFileSystemManager.buildConfigMap(exampleSas, container)
  val exampleStorageEndpoint = BlobPathBuilderSpec.buildEndpoint("testStorageAccount")
  val exampleCombinedEndpoint = BlobFileSystemManager.combinedEnpointContainerUri(exampleStorageEndpoint, container)


  // This test is passing locally but not in CI, because it can't find a filesystem provider with the azb scheme.
  // Some sort of build issue?
  it should "parse an expiration from a sas token" in {
    val provider = new AzureFileSystemProvider()
    val filesystem : AzureFileSystem = provider.newFileSystem(exampleCombinedEndpoint, exampleConfig.asJava).asInstanceOf[AzureFileSystem]
    filesystem.getExpiry.asScala shouldBe Some(now)
    filesystem.getFileStores.asScala.map(_.name()).exists(_ == container.value) shouldBe true
  }
}
