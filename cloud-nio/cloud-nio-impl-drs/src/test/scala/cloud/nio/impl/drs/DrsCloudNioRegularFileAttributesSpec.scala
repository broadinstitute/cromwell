package cloud.nio.impl.drs

import com.google.cloud.NoCredentials
import org.apache.http.impl.client.HttpClientBuilder
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration.Duration

class DrsCloudNioRegularFileAttributesSpec extends FlatSpecLike with Matchers {

  private val fakeDrsConfig = DrsConfig("http://martha-url", "{}")
  private val fakeCredentials = NoCredentials.getInstance
  private val httpClientBuilder = HttpClientBuilder.create()
  private val accessTokenAcceptableTTL: Duration = Duration.Inf

//  private def drsReadInterpreter(marthaResponse: MarthaResponse): IO[ReadableByteChannel] =
//    throw new UnsupportedOperationException("Currently DrsResolverSpec doesn't need to use drs read interpreter.")

  private val mockDrsPathResolver = new MockEngineDrsPathResolver(fakeDrsConfig, httpClientBuilder, fakeCredentials, accessTokenAcceptableTTL)
  private val mockDrsFileAttributesWithAllHashes = new MockDrsCloudNioRegularFileAttributes(MockDrsPaths.drsPathResolvingGcsPathWithAllHashes, mockDrsPathResolver)
  private val mockDrsFileAttributesWithMd5Hash = new MockDrsCloudNioRegularFileAttributes(MockDrsPaths.drsPathResolvingGcsPathWithMd5Hash, mockDrsPathResolver)
  private val mockDrsFileAttributesWithShaHash = new MockDrsCloudNioRegularFileAttributes(MockDrsPaths.drsPathResolvingGcsPathWithShaHash, mockDrsPathResolver)
  private val mockDrsFileAttributesWithNoPreferredHash = new MockDrsCloudNioRegularFileAttributes(MockDrsPaths.drsPathResolvingGcsPathWithNoPreferredHash, mockDrsPathResolver)

  behavior of "fileHash()"

  it should "return crc32c hash from `hashes` in Martha response when there is a crc32c" in {
    mockDrsFileAttributesWithAllHashes.fileHash shouldBe Option("8a366443")
  }

  it should "return md5 hash from `hashes` in Martha response when there is no crc32c" in {
    mockDrsFileAttributesWithMd5Hash.fileHash shouldBe Option("336ea55913bc261b72875bd259753046")
  }

  it should "return sha256 hash from `hashes` in Martha response when there is only a sha256" in {
    mockDrsFileAttributesWithShaHash.fileHash shouldBe Option("f76877f8e86ec3932fd2ae04239fbabb8c90199dab0019ae55fa42b31c314c44")
  }

  it should "return first (alphabetized by type) hash from `hashes` in Martha response when there are no preferred hash types" in {
    mockDrsFileAttributesWithNoPreferredHash.fileHash shouldBe Option("xrd")
  }
}
