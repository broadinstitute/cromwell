package cromwell.util.google

import com.google.api.services.storage.StorageScopes

object GoogleScopes {
  val Scopes = Vector( // FIXME: Should be in package object? I believe so
    "https://www.googleapis.com/auth/genomics",
    "https://www.googleapis.com/auth/devstorage.full_control",
    "https://www.googleapis.com/auth/devstorage.read_write",
    "https://www.googleapis.com/auth/compute",
    StorageScopes.DEVSTORAGE_READ_WRITE
  )
}
