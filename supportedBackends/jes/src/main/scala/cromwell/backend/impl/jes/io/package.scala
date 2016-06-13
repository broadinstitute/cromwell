package cromwell.backend.impl.jes

import java.nio.file.{Files, Path}

import com.google.api.client.http.HttpResponseException
import cromwell.backend.BackendWorkflowDescriptor
import cromwell.backend.impl.jes.JesImplicits.GoogleAuthWorkflowOptions
import cromwell.filesystems.gcs._

package object io {
  implicit class PathEnhanced(val path: Path) extends AnyVal {
    import better.files._

    def hash = path match {
      case gcs: NioGcsPath => gcs.getFileSystem.provider().asInstanceOf[GcsFileSystemProvider].crc32cHash(gcs)
      case _ => path.md5
    }

    def writeAsJson(content: String): File = {
      Files.write(path, content.getBytes, ContentTypeOption.Json)
    }
  }

  private [jes] def buildFilesystem(workflowDescriptor: BackendWorkflowDescriptor, authMode: GoogleAuthMode, googleConfig: GoogleConfiguration) = {
    val authOptions = workflowDescriptor.workflowOptions.toGoogleAuthOptions
    val storage = authMode.buildStorage(authOptions, googleConfig)

    GcsFileSystemProvider(storage).getFileSystem
  }

  private [jes] def isFatalJesException(t: Throwable): Boolean = t match {
    case e: HttpResponseException if e.getStatusCode == 403 => true
    case e: HttpResponseException if e.getStatusCode == 400 && e.getContent.contains("INVALID_ARGUMENT") => true
    case _ => false
  }

  private [jes] def isTransientJesException(t: Throwable): Boolean = t match {
    // Quota exceeded
    case e: HttpResponseException if e.getStatusCode == 429 => true
    case _ => false
  }
}
