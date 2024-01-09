package cromwell

import wdl.draft2.model._
import wom.JobOutput
import wom.values.WomValue

package object engine {

  implicit class EnhancedFullyQualifiedName(val fqn: FullyQualifiedName) extends AnyVal {
    def scopeAndVariableName: (String, String) = {
      val array = fqn.split("\\.(?=[^\\.]+$)")
      (array(0), array(1))
    }
  }

  implicit class EnhancedCallOutputMap[A](val m: Map[A, JobOutput]) extends AnyVal {
    def mapToValues: Map[A, WomValue] = m map { case (k, JobOutput(womValue)) =>
      (k, womValue)
    }
  }

}
