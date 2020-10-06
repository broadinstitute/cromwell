package cloud.nio.impl.drs

import cats.effect.IO
import cloud.nio.impl.drs.DrsCloudNioFileProvider.DrsReadInterpreter
import cloud.nio.impl.drs.MockDrsCloudNioFileSystemProvider._
import com.google.cloud.NoCredentials
import com.typesafe.config.{Config, ConfigFactory}
import org.apache.http.impl.client.HttpClientBuilder
import org.specs2.mock.Mockito
import org.specs2.mock.Mockito._

import scala.concurrent.duration.Duration

class MockDrsCloudNioFileSystemProvider(config: Config = mockConfig,
                                        httpClientBuilder: HttpClientBuilder = Mockito.mock[HttpClientBuilder].smart,
                                        drsReadInterpreter: DrsReadInterpreter = (_, _) =>
                                          IO.raiseError(
                                            new UnsupportedOperationException("mock did not specify a read interpreter")
                                          ),
                                        mockResolver: Option[EngineDrsPathResolver] = None,
                                       )
  extends DrsCloudNioFileSystemProvider(config, NoCredentials.getInstance, httpClientBuilder, drsReadInterpreter) {

  override lazy val drsPathResolver: EngineDrsPathResolver = {
    mockResolver getOrElse
      new MockEngineDrsPathResolver(
        drsConfig = drsConfig,
        httpClientBuilder = httpClientBuilder,
        accessTokenAcceptableTTL = Duration.Inf,
      )
  }
}

object MockDrsCloudNioFileSystemProvider {
  private lazy val mockConfig = ConfigFactory.parseString(
    """martha.url = "https://mock.martha"
      |access-token-acceptable-ttl = 1 hour
      |""".stripMargin
  )
}
