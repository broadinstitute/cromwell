package cromwell.engine.backend.jes

import cromwell.engine.backend.jes.authentication.JesDockerCredentials
import cromwell.engine.io.gcs.GoogleConfiguration
import cromwell.util.DockerConfiguration

trait JesConfiguration {
  def jesConf: JesAttributes
  def googleConf: GoogleConfiguration
  def dockerConf: Option[JesDockerCredentials]
}

object ProductionJesConfiguration {
  lazy val jesConf =  JesAttributes()
  lazy val googleConf =  GoogleConfiguration.gcloudConf.get
  lazy val dockerConf = DockerConfiguration.dockerConf.dockerCredentials map JesDockerCredentials.apply
}

trait ProductionJesConfiguration extends JesConfiguration {
  override lazy val jesConf = ProductionJesConfiguration.jesConf
  override lazy val googleConf = ProductionJesConfiguration.googleConf
  override lazy val dockerConf = ProductionJesConfiguration.dockerConf
}