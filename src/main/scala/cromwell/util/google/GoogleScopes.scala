package cromwell.util.google

import com.google.api.services.genomics.GenomicsScopes
import com.google.api.services.storage.StorageScopes
import scala.collection.JavaConverters._

object GoogleScopes {
  val Scopes = {
    val scopes = Vector(
      StorageScopes.DEVSTORAGE_FULL_CONTROL,
      StorageScopes.DEVSTORAGE_READ_WRITE,
      "https://www.googleapis.com/auth/genomics",
      "https://www.googleapis.com/auth/compute"
    ) ++ GenomicsScopes.all().asScala
    scopes.sorted.distinct
  }
}
