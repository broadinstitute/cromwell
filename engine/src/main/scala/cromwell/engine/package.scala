package cromwell

import cromwell.core.JobOutput
import wdl4s._
import wdl4s.values.WdlValue

package object engine {

  implicit class EnhancedFullyQualifiedName(val fqn: FullyQualifiedName) extends AnyVal {
    def scopeAndVariableName: (String, String) = {
      val array = fqn.split("\\.(?=[^\\.]+$)")
      (array(0), array(1))
    }
  }

  implicit class EnhancedCallOutputMap[A](val m: Map[A, JobOutput]) extends AnyVal {
    def mapToValues: Map[A, WdlValue] = m map {
      case (k, JobOutput(wdlValue)) => (k, wdlValue)
    }
  }

}
