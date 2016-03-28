package cromwell.engine.backend.jes

import cromwell.engine.backend.jes.authentication.JesDockerCredentials
import cromwell.filesystems.gcs.{GoogleConfigurationAdapter, GoogleConfiguration}
import cromwell.filesystems.gcs.GoogleConfigurationAdapter
import cromwell.util.DockerConfiguration

trait JesConfiguration {
  def jesConf: JesAttributes
  def genomicsConf: GoogleConfiguration
  def gcsConf: GoogleConfiguration
  def dockerConf: Option[JesDockerCredentials]
}

object ProductionJesConfiguration {
  lazy val jesConf =  JesAttributes()
  private lazy val oldGoogleConf =  GoogleConfigurationAdapter.gcloudConf.get
  lazy val gcsConf = GoogleConfiguration(oldGoogleConf.appName, oldGoogleConf.userAuthMode.getOrElse(oldGoogleConf.cromwellAuthMode))
  lazy val genomicsConf = GoogleConfiguration(oldGoogleConf.appName, oldGoogleConf.cromwellAuthMode)
  lazy val dockerConf = DockerConfiguration.dockerConf.dockerCredentials map JesDockerCredentials.apply
}

trait ProductionJesConfiguration extends JesConfiguration {
  override lazy val jesConf = ProductionJesConfiguration.jesConf
  override lazy val genomicsConf = ProductionJesConfiguration.genomicsConf
  override lazy val gcsConf = ProductionJesConfiguration.gcsConf
  override lazy val dockerConf = ProductionJesConfiguration.dockerConf
}