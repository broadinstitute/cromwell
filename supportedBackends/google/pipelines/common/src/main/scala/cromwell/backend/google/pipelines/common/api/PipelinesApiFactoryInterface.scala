package cromwell.backend.google.pipelines.common.api

import com.google.auth.Credentials

abstract class PipelinesApiFactoryInterface {
  def fromCredentials(credentials: Credentials): PipelinesApiRequestFactory
}
