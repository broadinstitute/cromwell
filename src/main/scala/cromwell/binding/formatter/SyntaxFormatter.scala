package cromwell.binding.formatter

import cromwell.binding.command.{ParameterCommandPart, StringCommandPart, Command}
import cromwell.binding.types.WdlType

import scala.collection.JavaConverters._

import cromwell.binding._
import cromwell.parser.WdlParser.{Terminal, AstList, AstNode, Ast}

trait SyntaxHighlighter {
  def keyword(s: String): String = s
  def name(s: String): String = s
  def section(s: String): String = s
  def wdlType(t: WdlType): String = t.toWdlString
  def variable(s: String): String = s
  def alias(s: String): String = s
  def command(s: String): String = s
  def function(s: String): String = s
}

object NullSyntaxHighlighter extends SyntaxHighlighter

object AnsiSyntaxHighlighter extends SyntaxHighlighter {
  def highlight(string: String, color: Int) = s"\033[38;5;${color}m${string}\033[0m"
  override def keyword(s: String): String = highlight(s, 214)
  override def name(s: String): String = highlight(s, 253)
  override def section(s: String): String = s
  override def wdlType(t: WdlType): String = highlight(t.toWdlString, 33)
  override def variable(s: String): String = highlight(s, 112)
  override def alias(s: String): String = s
  override def command(s: String): String = s
  override def function(s: String): String = highlight(s, 13)
}

object HtmlSyntaxHighlighter extends SyntaxHighlighter {
  def wrap(s: String, cls: String) = s"""<span class="$cls">$s</span>"""
  override def keyword(s: String): String = wrap(s, "keyword")
  override def name(s: String): String = wrap(s, "name")
  override def section(s: String): String = wrap(s, "section")
  override def wdlType(t: WdlType): String = wrap(t.toWdlString, "type")
  override def variable(s: String): String = wrap(s, "variable")
  override def alias(s: String): String = wrap(s, "alias")
  override def command(s: String): String = wrap(s, "command")
  override def function(s: String): String = wrap(s, "function")
}

class SyntaxFormatter(highlighter: SyntaxHighlighter = NullSyntaxHighlighter) {
  val indent = 2
  def format(binding: WdlBinding): String = {
    val imports = binding.imports.map(formatImport) match {
      case v if v.size > 0 => v.mkString("\n") + "\n\n"
      case v => ""
    }
    val definitions = for(node <- binding.ast.getAttribute("definitions").asInstanceOf[AstList].asScala.toVector) yield {
      node match {
        case a:Ast if a.getName == "Workflow" => formatWorkflow(binding.workflows.head)
        case a:Ast if a.getName == "Task" => formatTask(binding.findTask(text(a.getAttribute("name"))).getOrElse {
          throw new UnsupportedOperationException("Shouldn't happen")
        })
      }
    }
    s"$imports${definitions.mkString("\n\n")}"
  }
  private def formatImport(imp: Import): String = {
    val namespace = imp.namespace.map{ns => s" as $ns"}.getOrElse("")
    s"${highlighter.keyword("import")} '${imp.uri}'$namespace"
  }
  private def formatTask(task: Task): String = {
    val outputs = if (task.outputs.nonEmpty) formatOutputs(task.outputs, 1) else ""
    val command = formatCommandSection(task.command, 1)
    val sections = List(command, outputs).filter(_.size > 0)
    val header = s"""${highlighter.keyword("task")} ${highlighter.name(task.name)} {
       |${sections.mkString("\n")}
       |}"""
     .stripMargin
    header
  }
  private def formatCommandSection(command: Command, level:Int): String = {
    val section = s"""${highlighter.section("command")} {
        |${formatCommand(command, level+1)}
        |}"""
    indentText(section.stripMargin, level)
  }
  private def formatCommand(command: Command, level:Int): String = {
    val abstractCommand = command.parts.map {
      case x:StringCommandPart => x
      case x:ParameterCommandPart => s"$${${highlighter.wdlType(x.wdlType)} ${highlighter.variable(x.name)}}"
    }
    indentText(highlighter.command(command.normalize(abstractCommand.mkString)), 1)
  }
  private def formatOutputs(outputs: Seq[TaskOutput], level:Int): String = {
    val section = s"""${highlighter.section("output")} {
        |${outputs.map(formatOutput(_, 1)).mkString("\n")}
        |}"""
    indentText(section.stripMargin, level)
  }
  private def formatOutput(output: TaskOutput, level:Int): String = {
    indentText(s"${highlighter.wdlType(output.wdlType)} ${highlighter.variable(output.name)} = ${output.expression.toString(highlighter)}", level)
  }
  private def formatWorkflow(workflow: Workflow): String = {
    s"""${highlighter.keyword("workflow")} ${highlighter.name(workflow.name)} {
        |${workflow.calls.map{formatCall(_, 1)}.mkString("\n")}
        |}""".stripMargin
  }
  private def formatCall(call: Call, level:Int): String = {
    val header = s"${highlighter.keyword("call")} ${highlighter.name(call.task.name)}${formatCallAlias(call)}"
    if (call.inputMappings.isEmpty) {
      indentText(header, level)
    } else {
      val inputString = call.inputMappings.map {case (k, v) =>
        s"$k=${v.toString(highlighter)}"
      }.mkString(", ")
      indentText(s"""$header {
         |  input: $inputString
         |}""".stripMargin, level)
    }
  }
  private def formatCallAlias(call: Call): String = {
    call.alias.map {a => s" as ${highlighter.alias(a)}"}.getOrElse("")
  }
  private def text(astNode: AstNode) = astNode.asInstanceOf[Terminal].getSourceString
  private def indentText(s: String, i: Int): String = {
    s.split("\n").map {" " * (i * indent) + _}.mkString("\n")
  }
}
