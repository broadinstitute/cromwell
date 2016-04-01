package cromwell.engine.backend.jes.authentication

import com.google.api.services.genomics.Genomics
import cromwell.core.WorkflowOptions
import cromwell.engine.backend.io.filesystem.gcs.{GcsFileSystem, GcsFileSystemProvider, StorageFactory}
import cromwell.engine.backend.jes._

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

  lazy val genomicsInterface = GenomicsFactory(googleConf.appName, jesConf.endpointUrl)
  lazy val cromwellGcsFileSystem = StorageFactory.cromwellAuthenticated map { cromwellStorage =>
    GcsFileSystemProvider(cromwellStorage).getFileSystem
  } getOrElse { throw new Exception("JES Backend requires a GCS configuration. No suitable configuration has been found.") }
}

trait ProductionJesAuthentication extends JesConnection {
  override lazy val genomicsInterface = ProductionJesConnection.genomicsInterface
  override lazy val cromwellGcsFileSystem = ProductionJesConnection.cromwellGcsFileSystem

  // User authenticated filesystem defaults back to cromwell authenticated if there is no user configuration
  override def userGcsFileSystem(options: WorkflowOptions) = {
    StorageFactory.userAuthenticated(options) map { storage => GcsFileSystemProvider(storage).getFileSystem } getOrElse cromwellGcsFileSystem
  }
}
