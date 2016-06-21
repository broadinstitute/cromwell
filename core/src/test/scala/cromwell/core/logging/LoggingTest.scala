package cromwell.core.logging

import java.util.regex.Pattern

object LoggingTest {
  def escapePattern(pattern: String) = Pattern.quote(pattern)
}
