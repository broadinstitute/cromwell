package cloud.nio.impl.drs

import com.google.cloud.NoCredentials
import org.apache.http.impl.client.HttpClientBuilder
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration.Duration

class DrsCloudNioRegularFileAttributesSpec extends FlatSpecLike with Matchers {

  val fakeDrsConfig = DrsConfig("http://martha-url", "{}")
  val fakeCredentials = NoCredentials.getInstance
  val httpClientBuilder = HttpClientBuilder.create()
  val accessTokenAcceptableTTL: Duration = Duration.Inf

//  private def drsReadInterpreter(marthaResponse: MarthaResponse): IO[ReadableByteChannel] =
//    throw new UnsupportedOperationException("Currently DrsResolverSpec doesn't need to use drs read interpreter.")

  val mockDrsPathResolver = new MockDrsPathResolver(fakeDrsConfig, httpClientBuilder, fakeCredentials, accessTokenAcceptableTTL)
  val mockDrsFileAttributes = new MockDrsCloudNioRegularFileAttributes(MockDrsPaths.drsPathResolvingGcsPath, mockDrsPathResolver)

  behavior of "fileHash()"

  it should "return proper hash from `hashes` in Martha response" in {
    mockDrsFileAttributes.fileHash shouldBe Option("8a366443")
  }
}
