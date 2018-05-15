package cromwell.backend

import shapeless.{:+:, CNil}
import wom.callable.AdHocValue

package object standard {
  object StandardAdHocValue {
    object AsAdHocValue {
      def unapply(arg: StandardAdHocValue): Option[AdHocValue] = arg.select[AdHocValue]
    }
    object AsLocalizedAdHocValue {
      def unapply(arg: StandardAdHocValue): Option[LocalizedAdHocValue] = arg.select[LocalizedAdHocValue]
    }
  }
  
  type StandardAdHocValue = AdHocValue :+: LocalizedAdHocValue :+: CNil
}
