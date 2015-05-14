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

object TerminalSyntaxHighlighter extends SyntaxHighlighter {
  override def keyword(s: String): String = s"\033[38;5;214m$s\033[0m"
  override def name(s: String): String = s"\033[38;5;253m$s\033[0m"
  override def section(s: String): String = s
  override def wdlType(t: WdlType): String = s"\033[38;5;33m${t.toWdlString}\033[0m"
  override def variable(s: String): String = s"\033[38;5;112m$s\033[0m"
  override def alias(s: String): String = s
  override def command(s: String): String = s //s"\033[48;5;243m$s\033[0m"
  override def function(s: String): String = s"\033[38;5;13m$s\033[0m"
}

object SyntaxFormatter {
  def highlighter = TerminalSyntaxHighlighter
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
    val namespace = imp.namespace match {
      case Some(ns) => s" as $ns"
      case None => ""
    }
    s"${highlighter.keyword("import")} '${imp.uri}'$namespace"
  }
  private def formatTask(task: Task): String = {
    val outputs = if(task.outputs.size > 0) formatOutputs(task.outputs, 1) else ""
    val command = formatCommandSection(task.command, 1)
    val sections = List(command, outputs)
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
  private def formatCommand(command: Command, level:Int) = {
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
    s"${in(level)}${highlighter.wdlType(output.wdlType)} ${highlighter.variable(output.name)} = ${output.expression.toString(highlighter)}"
  }
  private def formatWorkflow(workflow: Workflow): String = {
    val section = s"""${highlighter.keyword("workflow")} ${highlighter.name(workflow.name)} {
        |${workflow.calls.map{formatCall(_, 1)}.mkString("\n")}
        |}"""
    section.stripMargin
  }
  private def formatCall(call: Call, level:Int): String = {
    val header = s"${highlighter.keyword("call")} ${highlighter.name(call.task.name)}${formatCallAlias(call)}"
    if (call.inputMappings.size == 0) {
      indentText(header, level)
    } else {
      val inputString = call.inputMappings.map{ case (k, v) =>
        s"$k=${v.toString(highlighter)}"
      }.mkString(", ")
      indentText(s"""$header {
         |  input: $inputString
         |}""".stripMargin, level)
    }
  }
  private def formatCallAlias(call: Call): String = {
    call.alias match {
      case Some(s:String) => s" as ${highlighter.alias(s)}"
      case None => ""
    }
  }
  private def text(astNode: AstNode) = astNode.asInstanceOf[Terminal].getSourceString
  private def indentText(s: String, i: Int): String = {
    s.split("\n").map{in(i) + _}.mkString("\n")
  }
  private def in(i: Int): String = " " * (i * indent)
}
