package cloud.nio.impl.drs

import cats.effect.IO
import cloud.nio.impl.drs.DrsCloudNioFileProvider.DrsReadInterpreter
import cloud.nio.impl.drs.MockDrsCloudNioFileSystemProvider.MockConfig
import com.google.auth.oauth2.OAuth2Credentials
import com.google.cloud.NoCredentials
import com.typesafe.config.{Config, ConfigFactory}
import org.apache.http.impl.client.HttpClientBuilder
import org.specs2.mock.Mockito
import org.specs2.mock.Mockito._

import scala.concurrent.duration.Duration

class MockDrsCloudNioFileSystemProvider(config: Config = MockConfig,
                                        credentials: OAuth2Credentials = NoCredentials.getInstance,
                                        httpClientBuilder: HttpClientBuilder = Mockito.mock[HttpClientBuilder].smart,
                                        drsReadInterpreter: DrsReadInterpreter = (_, _) =>
                                          IO.raiseError(
                                            new UnsupportedOperationException("mock did not specify a read interpreter")
                                          ),
                                        mockResolver: Option[EngineDrsPathResolver] = None,
                                       )
  extends DrsCloudNioFileSystemProvider(config, credentials, httpClientBuilder, drsReadInterpreter) {

  override lazy val drsPathResolver: EngineDrsPathResolver = {
    mockResolver getOrElse new MockEngineDrsPathResolver(drsConfig, credentials, httpClientBuilder, Duration.Inf)
  }
}

object MockDrsCloudNioFileSystemProvider {
  private lazy val MockConfig = ConfigFactory.parseString(
    """martha.url = "https://mock.martha"
      |access-token-acceptable-ttl = 1 hour
      |""".stripMargin
  )
}
