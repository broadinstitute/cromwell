package cwl

import cats.syntax.traverse._
import cats.syntax.validated._
import cats.instances.list._
import common.validation.ErrorOr._
import wom.types.WomNothingType
import wom.values._

import scala.collection.JavaConverters._

/**
  * Partial copy-port of cwltool's expression.py.
  *
  * Minimizes refactor so that as expression.py continues to update the updates may be manually copied here also.
  *
  * Current omissions in `def interpolate`:
  * - `def evaluator`
  *   - Utility not copy-ported.
  *   - `def interpolate` instead ignores absence of InlineJavascriptRequirement and instead always uses fullJS=True
  * - `leaf = json.dumps(e, sort_keys=True)`
  *   - When the interpolated string is not just-a-single-expression there is undefined behavior of rendering
  *     non-WomPrimitive values into json.
  *   - Even if it did, does not sort the keys.
  *   - Thus we do not currently return character-for-character equality for these expressions with cwltool, however
  *     these cases are not explicitly defined in the CWL spec as far as I know.
  *
  * @see https://github.com/common-workflow-language/cwltool/blob/353dbed/cwltool/expression.py
  */
//noinspection ZeroIndexToHead,RemoveRedundantReturn
object ExpressionInterpolator {

  class SubstitutionException(message: String) extends RuntimeException(message)

  /**
    * Copy-port of expression.py's scanner.
    */
  private def scanner(scan: String): List[Int] = {
    val DEFAULT = 0
    val DOLLAR = 1
    val PAREN = 2
    val BRACE = 3
    val SINGLE_QUOTE = 4
    val DOUBLE_QUOTE = 5
    val BACKSLASH = 6

    var i = 0
    val stack = new java.util.Stack[Int]
    stack.push(DEFAULT)
    var start = 0
    while (i < scan.length) {
      val state = stack.peek
      val c = scan(i)

      if (state == DEFAULT) {
        if (c == '$') {
          stack.push(DOLLAR)
        } else if (c == '\\') {
          stack.push(BACKSLASH)
        }
      } else if (state == BACKSLASH) {
        stack.pop()
        if (stack.peek == DEFAULT) {
          return List(i - 1, i + 1)
        }
      } else if (state == DOLLAR) {
        if (c == '(') {
          start = i - 1
          stack.push(PAREN)
        } else if (c == '{') {
          start = i - 1
          stack.push(BRACE)
        } else {
          stack.pop()
        }
      } else if (state == PAREN) {
        if (c == '(') {
          stack.push(PAREN)
        } else if (c == ')') {
          stack.pop()
          if (stack.peek == DOLLAR) {
            return List(start, i + 1)
          }
        } else if (c == '\'') {
          stack.push(SINGLE_QUOTE)
        } else if (c == '"') {
          stack.push(DOUBLE_QUOTE)
        }
      } else if (state == BRACE) {
        if (c == '{') {
          stack.push(BRACE)
        } else if (c == '}') {
          stack.pop()
          if (stack.peek == DOLLAR) {
            return List(start, i + 1)
          }
        } else if (c == '\'') {
          stack.push(SINGLE_QUOTE)
        } else if (c == '"') {
          stack.push(DOUBLE_QUOTE)
        }
      } else if (state == SINGLE_QUOTE) {
        if (c == '\'') {
          stack.pop()
        } else if (c == '\\') {
          stack.push(BACKSLASH)
        }
      } else if (state == DOUBLE_QUOTE) {
        if (c == '"') {
          stack.pop()
        } else if (c == '\\') {
          stack.push(BACKSLASH)
        }
      }
      i += 1
    }
    if (stack.size > 1) {
      throw new SubstitutionException(
        "Substitution error, unfinished block starting at position %s: %s".format(start, scan.drop(start)))
    } else {
      return null
    }
  }

  /**
    * Copy-port of expression.py's interpolate.
    */
  def interpolate(scanParam: String,
                  evaluator: String => ErrorOr[WomValue],
                  strip_whitespace: Boolean = true): ErrorOr[WomValue] = {
    var scan = scanParam
    if (strip_whitespace) {
      scan = scan.trim
    }
    val parts = new java.util.Stack[ErrorOr[WomValue]]
    var w = scanner(scan)
    while (w != null) {
      parts.push(WomString(scan.slice(0, w(0))).valid)

      if (scan(w(0)) == '$') {
        val e = evaluator(scan.slice(w(0) + 1, w(1)))
        if (w(0) == 0 && w(1) == scan.length && parts.size <= 1) {
          return e
        }
        /*
        Original-ish:
          var leaf = json.dumps(e, sort_keys=True)
          if (leaf(0) == '"') {
            leaf = leaf.drop(1).dropRight(1)
          }
        Instead, push the raw WomValue. WomString's won't have been quoted as python's json.dumps does.
        Other WomPrimitive's should be fine. CWL using other structures (Maps, Arrays, etc.) may have problems.
         */
        val leaf = e
        parts.push(leaf)
      } else if (scan(w(0)) == '\\') {
        val e = scan(w(1) - 1)
        parts.push(WomString(e.toString).valid)
      }

      scan = scan.drop(w(1))
      w = scanner(scan)
    }
    parts.push(WomString(scan).valid)
    return parts.asScala.toList.sequence[ErrorOr, WomValue] map { list =>
      WomString(list.map{
        //we represent nulls as this type because Wom doesn't have a "null" value, but it does have a nothing type
        case WomOptionalValue(WomNothingType, None) => "null"
        case other => other.valueString
      }.mkString(""))
    }
  }
}
