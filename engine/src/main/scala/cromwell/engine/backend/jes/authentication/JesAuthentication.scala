package cromwell.engine.backend.jes.authentication

import com.google.api.services.genomics.Genomics
import cromwell.core.WorkflowOptions
import cromwell.engine.backend.jes._
import cromwell.filesystems.gcs.{GoogleConfiguration, GcsFileSystemProvider}
import cromwell.filesystems.gcs.{GcsFileSystem, StorageFactory}

import scala.language.postfixOps

/**
 * Trait for JesConnection
 */
trait JesConnection {
  def genomicsInterface: Genomics
  def userGcsFileSystem(options: WorkflowOptions): GcsFileSystem
  def cromwellGcsFileSystem: GcsFileSystem
}

object ProductionJesConnection {
  import ProductionJesConfiguration._

  lazy val genomicsInterface = GenomicsFactory(genomicsConf, jesConf.endpointUrl)
  lazy val cromwellGcsFileSystem = GcsFileSystemProvider(StorageFactory(genomicsConf)).getFileSystem
}

trait ProductionJesAuthentication extends JesConnection {
  import ProductionJesConfiguration._

  override lazy val genomicsInterface = ProductionJesConnection.genomicsInterface
  override lazy val cromwellGcsFileSystem = ProductionJesConnection.cromwellGcsFileSystem

  override def userGcsFileSystem(options: WorkflowOptions) = {
    val refreshToken = options.get(GoogleConfiguration.RefreshTokenOptionKey).toOption
    GcsFileSystemProvider(StorageFactory(gcsConf, refreshToken)).getFileSystem
  }
}
