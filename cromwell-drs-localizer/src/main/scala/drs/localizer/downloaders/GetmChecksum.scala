package drs.localizer.downloaders

import cloud.nio.impl.drs.AccessUrl
import drs.localizer.downloaders.AccessUrlDownloader.Hashes

sealed trait GetmChecksum {
  def algorithm: String
  def value: String
  def args: String = s"--checksum-algorithm '$algorithm' --checksum '$value'"
}
case class Md5(override val value: String) extends GetmChecksum {
  override def algorithm: String = "md5"
}
case class Crc32c(override val value: String) extends GetmChecksum {
  override def algorithm: String = "gs_crc32c"
}
case class AwsEtag(override val value: String) extends GetmChecksum {
  override def algorithm: String = "s3_etag"
}
case object Null extends GetmChecksum {
  override def algorithm: String = "null"
  override def value: String = "null"
}
// The `value` for `Unsupported` will be the named algorithm keys
case class Unsupported(override val value: String) extends GetmChecksum {
  override def algorithm: String = "null"
}

object GetmChecksum {
  def apply(hashes: Hashes, accessUrl: AccessUrl): GetmChecksum = {
    hashes match {
      case Some(hashes) if hashes.nonEmpty =>
        // `hashes` uses the Martha keys for these hash algorithms, which in turn are forwarded DRS providers' keys for
        // the algorithms. `getm` has its own notions of what these algorithms are called. For the specific case of
        // `md5` the algorithm names are the same between DRS providers and `getm`, but all of the other algorithm names
        // currently differ between DRS providers and `getm`.
        if (hashes.contains("md5")) {
          Md5(hashes("md5"))
        }
        else if (hashes.contains("crc32c")) {
          Crc32c(hashes("crc32c"))
        }
        // etags could be anything; only ask `getm` to check s3 etags if this actually looks like an s3 signed url.
        else if (hashes.contains("etag") && accessUrl.url.matches("^https://[^/]+\\.s3\\.amazonaws.com/.*")) {
          AwsEtag(hashes("etag"))
        }
        // Not pictured: sha256 which is observed at least in staging data, e.g. open access Kids First object
        // drs://dg.F82A1A:5b92382f-51d2-424d-9def-c9ac0ed8b807. At the time of this writing `getm` does not support
        // sha256 checksums. This is not expected to be an issue since DRS appears to be standardizing on md5 for
        // checksums and md5 checksums have so far appeared alongside sha256 when sha256 is present.
        else {
          // If this code were running in Cromwell this condition would probably merit a warning but the localizer
          // runs on the VM and at best can only complain to stderr. The `getm` algorithm of `null` is specified which
          // means "do not validate checksums" with the stringified contents of the hashes map as a value.
          Unsupported(value = hashes.keys.mkString(", "))
        }
      case _ => Null
    }
  }
}
