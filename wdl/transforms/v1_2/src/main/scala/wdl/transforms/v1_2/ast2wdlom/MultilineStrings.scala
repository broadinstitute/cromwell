package wdl.transforms.v1_2.ast2wdlom

import wdl.v1_2.parser.WdlParser.Ast
import wdl.transforms.base.ast2wdlom.{GenericAst, GenericAstList, GenericAstNode, GenericTerminal}

/**
 * WDL 1.2 multi-line strings (`<<< ... >>>`).
 *
 * The leading whitespace rules from the specification are applied here, and the result is presented to the shared
 * expression converter as an ordinary `StringLiteral` made of string, escape and placeholder pieces.
 *
 * https://github.com/openwdl/wdl/blob/v1.2.1/SPEC.md#multi-line-strings
 */
object MultilineStrings {

  sealed private trait Item
  final private case class Character(c: Char) extends Item
  // An escape sequence or an expression placeholder. Never counts as whitespace.
  final private case class Opaque(node: GenericAstNode) extends Item
  private case object LineContinuation extends Item

  private val Newline = Character('\n')

  def toStringLiteral(ast: Ast): GenericAst = {
    val pieces: Seq[GenericAstNode] = Option(ast.getAttribute("pieces")).map(V1_2GenericAstNode.apply) match {
      case Some(list: GenericAstList) => list.astNodeList
      case Some(single) => Seq(single)
      case None => Seq.empty
    }

    val items: Vector[Item] = pieces.toVector flatMap {
      case t: GenericTerminal if t.getTerminalStr == "string" =>
        t.getSourceString.replace("\r\n", "\n").toVector.map(Character.apply)
      case t: GenericTerminal if t.getTerminalStr == "line_continuation" => Vector(LineContinuation)
      case other => Vector(Opaque(other))
    }

    val position = pieces.headOption.flatMap(_.firstTerminal)
    val line = position.map(_.getLine).getOrElse(0)
    val column = position.map(_.getColumn).getOrElse(0)

    SyntheticGenericAst("StringLiteral",
                        Map("pieces" -> SyntheticGenericAstList(regroup(stripWhitespace(items), line, column)))
    )
  }

  // In the context of multi-line strings, whitespace refers to space and tab characters only
  private def isWhitespace(item: Item): Boolean = item match {
    case Character(' ') | Character('\t') => true
    case _ => false
  }

  private def stripWhitespace(items: Vector[Item]): Vector[Item] = {
    // 1. Remove all line continuations and the whitespace that follows them
    val continued = {
      val builder = Vector.newBuilder[Item]
      var i = 0
      while (i < items.length)
        items(i) match {
          case LineContinuation =>
            i += 1
            while (i < items.length && isWhitespace(items(i))) i += 1
          case other =>
            builder += other
            i += 1
        }
      builder.result()
    }

    // 2. Remove all whitespace following the opening `<<<`, up to and including a newline (if any)
    val leadingRemoved = {
      val dropped = continued.dropWhile(isWhitespace)
      if (dropped.headOption.contains(Newline)) dropped.tail else dropped
    }

    // 3. Remove all whitespace preceding the closing `>>>`, up to and including a newline (if any)
    val trailingRemoved = {
      val dropped = leadingRemoved.reverse.dropWhile(isWhitespace)
      (if (dropped.headOption.contains(Newline)) dropped.tail else dropped).reverse
    }

    // 4. Use the non-blank lines to determine the common leading whitespace
    val lines = splitLines(trailingRemoved)
    val nonBlankLines = lines.filter(_.exists(item => !isWhitespace(item)))
    val commonLeadingWhitespace =
      if (nonBlankLines.isEmpty) 0 else nonBlankLines.map(_.takeWhile(isWhitespace).length).min

    // 5. Remove the common leading whitespace from each line (including blank lines that contain whitespace)
    val stripped = lines.map { line =>
      line.drop(math.min(commonLeadingWhitespace, line.takeWhile(isWhitespace).length))
    }
    stripped.zipWithIndex flatMap { case (line, index) => if (index == 0) line else Newline +: line }
  }

  private def splitLines(items: Vector[Item]): Vector[Vector[Item]] =
    items.foldLeft(Vector(Vector.empty[Item])) {
      case (lines, Newline) => lines :+ Vector.empty
      case (lines, item) => lines.init :+ (lines.last :+ item)
    }

  private def regroup(items: Vector[Item], line: Int, column: Int): Vector[GenericAstNode] = {
    val builder = Vector.newBuilder[GenericAstNode]
    val text = new StringBuilder

    def flushText(): Unit =
      if (text.nonEmpty) {
        builder += SyntheticGenericTerminal("string", text.toString, line, column)
        text.clear()
      }

    items foreach {
      case Character(c) => text += c
      case Opaque(node) =>
        flushText()
        builder += node
      case LineContinuation => // Already removed
    }
    flushText()
    builder.result()
  }
}
