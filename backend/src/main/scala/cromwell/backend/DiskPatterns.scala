package cromwell.backend

import scala.util.matching.Regex

/**
  * These patterns define the two disk patterns that should be supported by all backends
  * (because they're used by the best practice workflows)
  */
object DiskPatterns {
  val Identifier = "[a-zA-Z0-9-_]+"
  val Directory = """/[^\s]+"""
  val Integer = "[1-9][0-9]*"

  // e.g. local-disk 20 SSD
  val WorkingDiskPattern: Regex = s"""local-disk\\s+($Integer)\\s+($Identifier)""".r

  // e.g. /some/mnt 20 SSD
  val MountedDiskPattern: Regex = s"""($Directory)\\s+($Integer)\\s+($Identifier)""".r
}
