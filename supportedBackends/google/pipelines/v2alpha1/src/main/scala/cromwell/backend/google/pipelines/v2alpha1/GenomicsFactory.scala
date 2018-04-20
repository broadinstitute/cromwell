package cromwell.backend.google.pipelines.v2alpha1

import java.net.URL

import com.google.auth.Credentials
import cromwell.backend.google.pipelines.common.api.{PipelinesApiFactoryInterface, PipelinesApiRequestFactory}
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode

case class GenomicsFactory(applicationName: String, authMode: GoogleAuthMode, endpointUrl: URL) extends PipelinesApiFactoryInterface {
  def fromCredentials(credentials: Credentials): PipelinesApiRequestFactory = ???
}
