package cromwell.backend.impl

import com.google.api.services.genomics.Genomics
import cromwell.backend.BackendInitializationData
import cromwell.filesystems.gcs.GcsFileSystem

package object jes {
  case class JesBackendInitializationData(backendFilesystem: GcsFileSystem, genomics: Genomics) extends BackendInitializationData
}
