package cromwell.backend.google.pipelines.common.api

import com.google.api.client.http.HttpRequestInitializer
import com.google.auth.Credentials
import com.google.auth.http.HttpCredentialsAdapter

/**
  * The interface provides a single method to build a PipelinesApiRequestFactory
  * There should be one PipelinesApiRequestFactory created per workflow.
  * That is because they need credentials and those can be different for each workflow
  * A PipelinesApiRequestFactory is able to generate all the API requests that are needed to run jobs on the
  * Pipelines API for this workflow.
  */
abstract class PipelinesApiFactoryInterface {
  final def fromCredentials(credentials: Credentials): PipelinesApiRequestFactory = {
    val httpCredentialsAdapter = new HttpCredentialsAdapter(credentials)
    build(httpCredentialsAdapter)
  }
  
  protected def build(httpRequestInitializer: HttpRequestInitializer): PipelinesApiRequestFactory

  def usesEncryptedDocker: Boolean
}
