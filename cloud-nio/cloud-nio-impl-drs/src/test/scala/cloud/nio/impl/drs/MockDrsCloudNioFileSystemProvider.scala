package cloud.nio.impl.drs

import cats.effect.IO
import cloud.nio.impl.drs.DrsCloudNioFileProvider.DrsReadInterpreter
import cloud.nio.impl.drs.MockDrsCloudNioFileSystemProvider._
import com.google.cloud.NoCredentials
import com.typesafe.config.{Config, ConfigFactory}
import org.apache.http.impl.client.HttpClientBuilder

class MockDrsCloudNioFileSystemProvider(config: Config = mockConfig,
                                        httpClientBuilder: Option[HttpClientBuilder] = None,
                                        drsReadInterpreter: DrsReadInterpreter = (_, _) =>
                                          IO.raiseError(
                                            new UnsupportedOperationException("mock did not specify a read interpreter")
                                          ),
                                        mockResolver: Option[DrsPathResolver] = None
) extends DrsCloudNioFileSystemProvider(config,
                                        GoogleOauthDrsCredentials(NoCredentials.getInstance, config),
                                        drsReadInterpreter
    ) {

  override lazy val drsPathResolver: DrsPathResolver =
    mockResolver getOrElse
      new MockDrsPathResolver(
        drsConfig = drsConfig,
        httpClientBuilderOverride = httpClientBuilder
      )
}

object MockDrsCloudNioFileSystemProvider {
  private lazy val mockConfig = ConfigFactory.parseString(
    """resolver.url = "https://mock.drshub"
      |access-token-acceptable-ttl = 1 hour
      |""".stripMargin
  )
}
