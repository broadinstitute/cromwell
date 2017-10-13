package cromwell

import wdl._
import wom.JobOutput
import wom.values.WdlValue

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
