package cromwell

import wdl4s.values.WdlValue

import scala.language.postfixOps
import scala.util.Success

package object backend {
  implicit class AugmentedAttemptedLookupSequence(s: Seq[AttemptedLookupResult]) {
    def toLookupMap: Map[String, WdlValue] = s collect {
      case AttemptedLookupResult(name, Success(value)) => (name, value)
    } toMap
  }
}
