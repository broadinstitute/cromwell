package cromwell.filesystems

import com.google.api.services.storage.StorageScopes

import scala.collection.JavaConverters._

package object gcs {
  type RefreshToken = String
  type GoogleScopes = java.util.Collection[String]

  val GcsScopes: GoogleScopes = List(
    StorageScopes.DEVSTORAGE_FULL_CONTROL,
    StorageScopes.DEVSTORAGE_READ_WRITE
  ).asJava
}
