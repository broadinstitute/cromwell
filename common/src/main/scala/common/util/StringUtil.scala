package common.util

import common.util.UriUtil._
import mouse.all._
import org.apache.commons.lang3.StringUtils
import org.apache.commons.text.StringEscapeUtils

import java.net.URI
import scala.util.Try

object StringUtil {

  implicit class EnhancedToStringable(val any: Any) extends AnyVal {
    // Uses pprint to (safely-ish) convert [Any] to a shortened string.
    // NB: the limit will be approximate due to the need to square-root it first (see below). You'll actually
    //   get the next highest square number as your limit.
    def toPrettyElidedString(limit: Int): String = {
      val sqrt = Math.max(Math.ceil(Math.sqrt(limit.doubleValue())).intValue(), 1)

      // I expected this to limit both width and height to the limit (eg a 100 wide by 100 high block of text).
      // In fact, it seems to limit total output length to a total of (width x height) characters.
      // ... shrug... ok fine. It's still a pretty output.
      pprint.apply(any, width = sqrt, height = sqrt).plainText
    }
  }

  implicit class EnhancedString(val string: String) extends AnyVal {

    /**
      * Ensure string ends with a /
      */
    def ensureSlashed: String = StringUtils.appendIfMissing(string, "/")

    /**
      * Ensure string does not end with a /
      */
    def ensureUnslashed: String = string.ensureSlashed |> StringUtils.chop

    /**
      * Ensure string does not start with a /
      */
    def ensureNoLeadingSlash: String = string.stripPrefix("/")

    /**
      * Escape for shell use
      */
    def ensureShellEscaped: String = StringEscapeUtils.escapeXSI(string)

    /**
      * Escape / with \/
      */
    def ensureSlashEscaped: String = string.replaceAll("/", "\\\\/")

    /**
      * Escape the string for use in shell and also escape slashes so it can be used in a sed expression while keeping
      * / as a sed separator
      */
    def ensureSedEscaped: String = string.ensureShellEscaped.ensureSlashEscaped

    /**
      * Makes the string look like a relative directory.
      * i.e no slash prefix and a slash suffix
      * e.g: /root/some/dir -> root/some/dir/
      */
    def relativeDirectory: String = string.ensureNoLeadingSlash.ensureSlashed

    def elided(limit: Int): String = {
      if (string.length > limit) {
        s"(elided) ${string.take(limit)}..."
      } else string
    }

    /**
      * Removes userInfo and sensitive query parts from strings that are RFC 2396 URIs.
      *
      * If the URI query does not contain any detected sensitive information, then the entire query will be masked.
      *
      * Returns the original input for:
      * - Any URI that cannot be parsed
      * - RFC 3986 URIs as those are not supported by java.net.URI
      * - Other non-IETF / non-W3C URI specifications such as DOS URI or DRS URI
      *
      * Depending on the encoding used in the input URI the masked output may have unexpected encoding. See:
      * - the StringUtilSpec for current expectations
      * - https://stackoverflow.com/questions/4571346/how-to-encode-url-to-avoid-special-characters-in-java#answer-4571518
      */
    def maskSensitiveUri: String = {
      Try(new URI(string))
        .map(_.maskSensitive.toASCIIString)
        .getOrElse(string)
    }
  }
}
