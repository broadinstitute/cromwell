package wdl.model.draft3.elements

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement.{AsciiCharacterEscape, BackslashEscape, DoubleQuoteEscape, NewlineEscape, SingleQuoteEscape, StringEscapeSequence, TabEscape, UnicodeCharacterEscape}

object StringEscapeSequence {
  val Octal = "\\\\([0-7]{3})".r
  val Hex = "\\\\x([0-9a-fA-F]{2})".r
  val FourDigitUnicode = "\\\\u([0-9a-fA-F]{4})".r
  val EightDigitUnicode = "\\\\U([0-9a-fA-F]{8})".r

  def parseEscapeSequence(seq: String): ErrorOr[StringEscapeSequence] = seq match {
    case "\\n" => NewlineEscape.validNel
    case "\\t" => TabEscape.validNel
    case "\\\"" => DoubleQuoteEscape.validNel
    case "\\'" => SingleQuoteEscape.validNel
    case "\\\\" => BackslashEscape.validNel
    case Octal(codePoint) => AsciiCharacterEscape(BigInt(codePoint, 8).byteValue()).validNel
    case Hex(codePoint) => AsciiCharacterEscape(BigInt(codePoint, 16).byteValue()).validNel
    case FourDigitUnicode(codePoint) => UnicodeCharacterEscape(BigInt(codePoint, 16).intValue).validNel
    case EightDigitUnicode(codePoint) => UnicodeCharacterEscape(BigInt(codePoint, 16).intValue).validNel


    case _ => s"Unrecognized escape sequence '$seq'".invalidNel
  }
}
