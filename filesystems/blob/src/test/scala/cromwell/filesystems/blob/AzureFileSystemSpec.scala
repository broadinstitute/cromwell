package cromwell.filesystems.blob

import com.azure.core.credential.AzureSasCredential
import com.azure.storage.blob.nio.{AzureFileSystem, AzureFileSystemProvider}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.time.{Duration, Instant}
import java.time.temporal.ChronoUnit
import scala.compat.java8.OptionConverters._
import scala.jdk.CollectionConverters._

class AzureFileSystemSpec extends AnyFlatSpec with Matchers {

  val fiveMinutes: Duration = Duration.of(5, ChronoUnit.MINUTES)

  private def makeFilesystemWithExpiration(expiration: Instant): AzureFileSystem =
    makeFilesystemWithCreds(BlobPathBuilderFactorySpec.buildExampleSasToken(expiration))

  private def makeFilesystemWithCreds(creds: AzureSasCredential): AzureFileSystem = {
    val storageEndpoint = BlobPathBuilderSpec.buildEndpoint("testStorageAccount")
    val container = BlobContainerName("testContainer")
    val combinedEndpoint = BlobFileSystemManager.combinedEnpointContainerUri(storageEndpoint, container)

    val provider = new AzureFileSystemProvider()
    provider
      .newFileSystem(
        combinedEndpoint,
        BlobFileSystemManager.buildConfigMap(creds, container).asJava
      )
      .asInstanceOf[AzureFileSystem]
  }

  it should "parse an expiration from a sas token" in {
    val now = Instant.now()
    val filesystem: AzureFileSystem = makeFilesystemWithExpiration(now)
    filesystem.getExpiry.asScala shouldBe Some(now)
    filesystem.getFileStores.asScala.map(_.name()).exists(_ == "testContainer") shouldBe true
  }

  it should "not be expired when the token is fresh" in {
    val anHourFromNow = Instant.now().plusSeconds(3600)
    val filesystem: AzureFileSystem = makeFilesystemWithExpiration(anHourFromNow)
    filesystem.isExpired(fiveMinutes) shouldBe false
  }

  it should "be expired when we're within the buffer" in {
    val threeMinutesFromNow = Instant.now().plusSeconds(180)
    val filesystem: AzureFileSystem = makeFilesystemWithExpiration(threeMinutesFromNow)
    filesystem.isExpired(fiveMinutes) shouldBe true
  }

  it should "be expired when the token is stale" in {
    val anHourAgo = Instant.now().minusSeconds(3600)
    val filesystem: AzureFileSystem = makeFilesystemWithExpiration(anHourAgo)
    filesystem.isExpired(fiveMinutes) shouldBe true
  }

  it should "not be expired with public credentials" in {
    val fileSystem = makeFilesystemWithCreds(BlobFileSystemManager.PLACEHOLDER_TOKEN)
    fileSystem.isExpired(fiveMinutes) shouldBe false
  }
}
