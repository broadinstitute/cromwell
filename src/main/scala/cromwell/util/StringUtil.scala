package cromwell.util

import org.apache.commons.codec.digest.DigestUtils

object StringUtil {

  implicit class EnhancedString(val value: String) extends AnyVal with Hashable {
    def md5Sum: String = DigestUtils.md5Hex(value)
  }

}
