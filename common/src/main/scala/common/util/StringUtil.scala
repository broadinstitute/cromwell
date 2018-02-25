package common.util

import mouse.all._
import org.apache.commons.lang3.StringUtils

object StringUtil {
  implicit class EnhancedString(val string: String) extends AnyVal {

    def ensureSlashed: String = StringUtils.appendIfMissing(string, "/")

    def ensureUnslashed: String = string.ensureSlashed |> StringUtils.chop
  }
}
