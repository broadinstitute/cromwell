package cromwell.util.google

import scala.collection.JavaConverters._

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object GenomicsScopes {
  val genomicsScopes = List(
    "https://www.googleapis.com/auth/genomics",
    "https://www.googleapis.com/auth/compute"
  ).asJava
}
