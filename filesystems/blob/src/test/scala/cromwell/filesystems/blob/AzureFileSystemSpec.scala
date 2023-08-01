package cromwell.filesystems.blob

import com.azure.storage.blob.nio.{AzureFileSystem, AzureFileSystemProvider}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.FileSystems
import java.nio.file.spi.FileSystemProvider
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
 // FileSystemProvider.installedProviders()

  it should "parse an expiration from a sas token" in {
    new AzureFileSystemProvider().getScheme shouldBe "azb"
    FileSystemProvider.installedProviders().asScala.map(_.getScheme) shouldBe List("azb")
    val fs = FileSystems.newFileSystem(exampleCombinedEndpoint, exampleConfig.asJava).asInstanceOf[AzureFileSystem]
    fs.getExpiry.asScala shouldBe Some(now)
    fs.getFileStores.asScala.map(_.name()).exists(_ == container.value) shouldBe true
  }
}
