package cromwell.filesystems

import scalaz.ValidationNel

package object gcs {
  type ErrorOr[+A] = ValidationNel[String, A]
  type RefreshToken = String
}
