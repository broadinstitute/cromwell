package cromwell.backend

import org.scalacheck.Properties
import org.scalacheck.Prop._

class RuntimeAttributeValidationSpec extends Properties("Runtime Validation") {

  property("foo") = secure{ true}


}
