package cromwell.util.google

import scala.collection.JavaConverters._

object GenomicsScopes {
  val genomicsScopes = List(
    "https://www.googleapis.com/auth/genomics",
    "https://www.googleapis.com/auth/compute"
  ).asJava
}
