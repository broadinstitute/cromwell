package cromwell.util.google

import com.google.api.services.storage.StorageScopes

object GoogleScopes {
  val Scopes = Vector(
    StorageScopes.DEVSTORAGE_FULL_CONTROL,
    StorageScopes.DEVSTORAGE_READ_WRITE,
    "https://www.googleapis.com/auth/genomics",
    "https://www.googleapis.com/auth/compute"
  )
}
