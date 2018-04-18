package wdl.draft3.transforms.ast2wdlom

import common.Checked
import wdl.draft3.parser.WdlParser.Ast
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.CommandPartElement.StringCommandPartElement
import wdl.model.draft3.elements.{CommandPartElement, CommandSectionElement, CommandSectionLine}

object AstToCommandSectionElement {
  def convert(ast: Ast): Checked[CommandSectionElement] = {
    ast.getAttributeAsVector[CommandPartElement]("parts") map { parts =>
      val lines = makeLines(parts)
      val trimmed = trimStartAndEnd(lines)
      val commonPrefixLength = minimumStartLength(trimmed)
      val commonPrefix = " " * commonPrefixLength

      CommandSectionElement(stripStarts(trimmed, commonPrefix))
    }
  }

  private def makeLines(elements: Vector[CommandPartElement]): Vector[CommandSectionLine] = {

    var accumulator = Vector.empty[CommandSectionLine]
    var current = Vector.empty[CommandPartElement]
    elements foreach {
      case StringCommandPartElement(str) if str.contains(System.lineSeparator) =>
        val split = str.split(System.lineSeparator, -1)
        accumulator :+= CommandSectionLine(current :+ StringCommandPartElement(split.head))
        accumulator ++= split.tail.init.map(s => CommandSectionLine(Vector(StringCommandPartElement(s))))
        current = split.tail.lastOption.map(StringCommandPartElement).toVector
      case other => current :+= other
    }
    val finalAccumulator = if (current.nonEmpty) {
      accumulator :+ CommandSectionLine(current)
    } else {
      accumulator
    }

    finalAccumulator map dropEmpties
  }

  private def trimStartAndEnd(elements: Vector[CommandSectionLine]): Vector[CommandSectionLine] = {
    elements.dropWhile(allWhitespace).reverse.dropWhile(allWhitespace).reverse
  }

  private def dropEmpties(line: CommandSectionLine): CommandSectionLine = {
    def empty(c: CommandPartElement): Boolean = c match {
      case StringCommandPartElement(s) if s.isEmpty => true
      case _ => false
    }
    CommandSectionLine(line.parts.filterNot(empty))
  }

  private def minimumStartLength(lines: Vector[CommandSectionLine]): Int = {
    val r = "^(\\s*)".r
    def startLength(line: CommandSectionLine) = line.parts.headOption match {
      case Some(StringCommandPartElement(str)) => r.findFirstIn(str).map(_.length).getOrElse(0)
      case _ => 0
    }

    lines.map(startLength).min
  }

  private def stripStarts(lines: Vector[CommandSectionLine], prefix: String): Vector[CommandSectionLine] = {
    if (prefix.isEmpty) lines else lines map { line =>
      line.parts.headOption match {
        case Some(StringCommandPartElement(str)) if str.startsWith(prefix) =>
          CommandSectionLine(Vector(StringCommandPartElement(str.stripPrefix(prefix))) ++ line.parts.tail)
        case _ => ??? // TODO: Should be impossible. Can we make the function total by fixing the inputs?
      }
    }
  }

  private def allWhitespace(s: String): Boolean = s.forall(_.isWhitespace)
  private def allWhitespace(c: CommandPartElement): Boolean = c match {
    case StringCommandPartElement(s) => allWhitespace(s)
    case _ => false
  }
  private def allWhitespace(l: CommandSectionLine): Boolean = l.parts.forall(allWhitespace)
}
