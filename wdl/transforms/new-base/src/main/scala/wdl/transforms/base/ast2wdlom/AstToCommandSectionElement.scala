package wdl.transforms.base.ast2wdlom

import cats.implicits._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.Checked._
import org.apache.commons.text.StringEscapeUtils
import wdl.model.draft3.elements.CommandPartElement.StringCommandPartElement
import wdl.model.draft3.elements.{CommandPartElement, CommandSectionElement, CommandSectionLine}

object AstToCommandSectionElement {
  def astToCommandSectionElement(implicit
    astNodeToCommandPartElement: CheckedAtoB[GenericAstNode, CommandPartElement]
  ): CheckedAtoB[GenericAst, CommandSectionElement] = CheckedAtoB.fromCheck { ast: GenericAst =>
    ast.getAttributeAsVector[CommandPartElement]("parts") flatMap { parts =>
      val lines = makeLines(parts)
      val trimmed = trimStartAndEndBlankLines(lines)

      // Commands support either tabs or spaces. For the undefined case of mixed tabs and spaces, Cromwell returns an error.
      // https://github.com/openwdl/wdl/blob/master/versions/1.0/SPEC.md#stripping-leading-whitespace
      val leadingWhitespaceMap = leadingWhitespace(trimmed)
      val distinctLeadingWhitespaceCharacters = leadingWhitespaceMap.mkString.distinct

      val commonPrefix: Checked[String] = distinctLeadingWhitespaceCharacters.length match {
        case 0 => "".validNelCheck
        case 1 =>
          (distinctLeadingWhitespaceCharacters.head.toString * leadingWhitespaceMap.map(_.length).min).validNelCheck
        case _ =>
          val charList = distinctLeadingWhitespaceCharacters.map { c: Char =>
            "\"" + StringEscapeUtils.escapeJava(c.toString) + "\""
          }

          s"Cannot mix leading whitespace characters in command: [${charList.mkString(", ")}]".invalidNelCheck
      }

      for {
        prefix <- commonPrefix
        lines <- stripStarts(trimmed, prefix)
      } yield CommandSectionElement(lines)
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

  private def trimStartAndEndBlankLines(elements: Vector[CommandSectionLine]): Vector[CommandSectionLine] =
    elements.dropWhile(allWhitespace).reverse.dropWhile(allWhitespace).reverse

  private def dropEmpties(line: CommandSectionLine): CommandSectionLine = {
    def empty(c: CommandPartElement): Boolean = c match {
      case StringCommandPartElement(s) if s.isEmpty => true
      case _ => false
    }
    CommandSectionLine(line.parts.filterNot(empty))
  }

  private def leadingWhitespace(lines: Vector[CommandSectionLine]): Vector[String] = {
    val r = "^(\\s*)".r

    def leadingWhitespaceForLine(line: CommandSectionLine): Option[String] = line.parts.headOption match {
      case Some(StringCommandPartElement(str)) => r.findFirstIn(str)
      case _ => None
    }

    lines.map(leadingWhitespaceForLine(_).getOrElse(""))
  }

  private def stripStarts(lines: Vector[CommandSectionLine], prefix: String): Checked[List[CommandSectionLine]] =
    if (prefix.isEmpty)
      lines.toList.validNelCheck
    else
      lines.toList traverse { line: CommandSectionLine =>
        line.parts.headOption match {
          case Some(StringCommandPartElement(str)) if str.startsWith(prefix) =>
            CommandSectionLine(
              Vector(StringCommandPartElement(str.stripPrefix(prefix))) ++ line.parts.tail
            ).validNelCheck
          case _ =>
            "Failed to strip common whitespace prefix from line.".invalidNelCheck
        }
      }

  private def allWhitespace(s: String): Boolean = s.forall(_.isWhitespace)
  private def allWhitespace(c: CommandPartElement): Boolean = c match {
    case StringCommandPartElement(s) => allWhitespace(s)
    case _ => false
  }
  private def allWhitespace(l: CommandSectionLine): Boolean = l.parts.forall(allWhitespace)
}
