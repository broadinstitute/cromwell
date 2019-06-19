package wdl.draft2.model.formatter

import common.util.TerminalUtil
import wdl.draft2.model.AstTools.EnhancedAstNode
import wdl.draft2.model._
import wdl.draft2.model.command.StringCommandPart
import wdl.draft2.parser.WdlParser.{Ast, AstList, AstNode}
import wom.types.WomType

import scala.collection.JavaConverters._

trait SyntaxHighlighter {
  def keyword(s: String): String = s
  def name(s: String): String = s
  def section(s: String): String = s
  def womType(t: WomType): String = t.stableName
  def variable(s: String): String = s
  def alias(s: String): String = s
  def command(s: String): String = s
  def function(s: String): String = s
}

object NullSyntaxHighlighter extends SyntaxHighlighter

object AnsiSyntaxHighlighter extends SyntaxHighlighter {
  override def keyword(s: String): String = TerminalUtil.highlight(214, s)
  override def name(s: String): String = TerminalUtil.highlight(253, s)
  override def section(s: String): String = keyword(s)
  override def womType(t: WomType): String = TerminalUtil.highlight(33, t.stableName)
  override def variable(s: String): String = TerminalUtil.highlight(112, s)
  override def alias(s: String): String = s
  override def command(s: String): String = s
  override def function(s: String): String = TerminalUtil.highlight(13, s)
}

object HtmlSyntaxHighlighter extends SyntaxHighlighter {
  def wrap(s: String, cls: String) = s"""<span class="$cls">$s</span>"""
  override def keyword(s: String): String = wrap(s, "keyword")
  override def name(s: String): String = wrap(s, "name")
  override def section(s: String): String = wrap(s, "section")
  override def womType(t: WomType): String = wrap(t.stableName, "type")
  override def variable(s: String): String = wrap(s, "variable")
  override def alias(s: String): String = wrap(s, "alias")
  override def command(s: String): String = wrap(s, "command")
  override def function(s: String): String = wrap(s, "function")
}

class SyntaxFormatter(highlighter: SyntaxHighlighter = NullSyntaxHighlighter) {
  val indentLevel = 2

  private def indent(s: String, i: Int): String = {
    s.split("\n").map {" " * (i * indentLevel) + _}.mkString("\n")
  }

  def format(namespace: WdlNamespace): String = {
    val imports = namespace.imports.map(formatImport) match {
      case v if v.nonEmpty => v.mkString("\n") + "\n\n"
      case _ => ""
    }

    /*
     TODO/FIXME: If 'body' is really a function of `namespace` then `WdlNamespace` should have a func which
     does the first part, and `NamespaceWithWorkflow` override it, call super and then adds on the second part
    */

    val namespaceDefinitions = namespace.ast.getAttribute("body").asInstanceOf[AstList].asScala.toVector

    val taskDefinitions = namespaceDefinitions collect { case a: Ast if a.getName == "Task" =>
      formatTask(namespace.findTask(a.getAttribute("name").sourceString).getOrElse(throw new UnsupportedOperationException("Shouldn't happen")))
    }

    val workflowDefinitions = namespace match {
      case n: WdlNamespaceWithWorkflow => namespaceDefinitions collect {case a: Ast if a.getName == "Workflow" => formatWorkflow(n.workflow)}
      case _ => Vector.empty[AstNode]
    }
    val definitions = taskDefinitions ++ workflowDefinitions

    s"$imports${definitions.mkString("\n\n")}"
  }

  private def formatImport(imp: Import): String = {
    val namespace = s" as ${imp.namespaceTerminal.sourceString}"

    s"${highlighter.keyword("import")} '${imp.uri}'$namespace"
  }

  private def formatTask(task: WdlTask): String = {
    val outputs = if (task.outputs.nonEmpty) formatOutputs(task.outputs, 1) else ""
    val command = formatCommandSection(task, 1)
    val declarations = task.declarations.map(formatDeclaration(_, 1)) match {
      case x: Seq[String] if x.nonEmpty => x.mkString("\n")
      case _ => ""
    }
    val runtime = formatRuntimeSection(task.runtimeAttributes, 1)
    val meta = formatMetaSection("meta", task.meta, 1)
    val parameterMeta = formatMetaSection("parameter_meta", task.parameterMeta, 1)
    val sections = List(declarations, command, outputs, runtime, meta, parameterMeta).filter(_.nonEmpty)
    val header = s"""${highlighter.keyword("task")} ${highlighter.name(task.name)} {
       |${sections.mkString("\n")}
       |}"""
     .stripMargin
    header
  }

  private def formatMetaSection(section: String, attrs: Map[String, String], level: Int): String = {
    attrs match {
      case m: Map[String, String] if m.nonEmpty =>
        val wdlAttrs = m map { case (k, v) => indent(s"$k: " + "\"" + v + "\"", 1) }
        indent(
          s"""${highlighter.keyword(section)} {
             |${wdlAttrs.mkString("\n")}
             |}""".stripMargin, level)
      case _ => ""
    }
  }

  private def formatRuntimeSection(runtimeAttributes: WdlRuntimeAttributes, level: Int): String = {
    runtimeAttributes.attrs match {
      case m if m.nonEmpty =>
        val attrs = m map { case (k, v) => indent(s"$k: ${v.toWomString}", level) }
        indent(
          s"""${highlighter.keyword("runtime")} {
             |${attrs.mkString("\n")}
             |}""".stripMargin, level)
      case _ => ""
    }
  }

  private def formatCommandSection(task: WdlTask, level:Int): String = {
    val (sdelim: String, edelim: String) =
      if (task.commandTemplate.collect({case s:StringCommandPart => s.literal}).mkString.contains("}")) ("<<<", ">>>")
      else ("{", "}")

    val section = s"""${highlighter.section("command")} $sdelim
        |${indent(highlighter.command(task.commandTemplateString), 1)}
        |$edelim"""
    indent(section.stripMargin, level)
  }

  private def formatOutputs(outputs: Seq[TaskOutput], level:Int): String = {
    val section = s"""${highlighter.section("output")} {
        |${outputs.map(formatOutput(_, 1)).mkString("\n")}
        |}"""
    indent(section.stripMargin, level)
  }

  private def formatOutput(output: TaskOutput, level:Int): String = {
    indent(s"${highlighter.womType(output.womType)} ${highlighter.variable(output.unqualifiedName)} = ${output.requiredExpression.toString(highlighter)}", level)
  }

  private def formatWorkflow(workflow: WdlWorkflow): String = {
    val declarations = workflow.declarations.map(formatDeclaration(_, 1))
    val children = workflow.children.collect({case c if !workflow.declarations.contains(c) => formatScope(c, 1) })
    val outputs = formatWorkflowOutputs(workflow.workflowOutputWildcards, 1)
    val meta = formatMetaSection("meta", workflow.meta, 1)
    val parameterMeta = formatMetaSection("parameter_meta", workflow.parameterMeta, 1)
    val sections = (declarations ++ children ++ Seq(meta, parameterMeta, outputs)).filter(_.nonEmpty)
    s"""${highlighter.keyword("workflow")} ${highlighter.name(workflow.unqualifiedName)} {
        |${sections.mkString("\n")}
        |}""".stripMargin
  }

  private def formatWorkflowOutputs(outputs: Seq[WorkflowOutputWildcard], level: Int): String = {
    outputs match {
      case x: Seq[WorkflowOutputWildcard] if x.nonEmpty =>
        val outputStrings = outputs.map(formatWorkflowOutput(_, 1))
        indent(s"""${highlighter.keyword("output")} {
                  |${outputStrings.mkString("\n")}
                  |}""".stripMargin, level)
      case _ => ""
    }
  }

  private def formatWorkflowOutput(output: WorkflowOutputWildcard, level: Int): String = {
    output.wildcard match {
      case true => indent(s"${formatWorkflowOutputFqn(output.fqn)}.*", level)
      case false => indent(formatWorkflowOutputFqn(output.fqn), level)
    }
  }

  private def formatWorkflowOutputFqn(fqn: String) = fqn.replaceFirst("[a-zA-Z0-9]+\\.", "")

  private def formatDeclaration(decl: DeclarationInterface, level: Int): String = {
    val expression = decl.expression.map(e => s" = ${e.toWomString}").getOrElse("")
    indent(s"${highlighter.womType(decl.womType)} ${highlighter.variable(decl.unqualifiedName)}$expression", level)
  }

  private def formatScope(scope: Scope, level: Int): String = scope match {
    case c: WdlTaskCall => formatCall(c, level)
    case s: Scatter => formatScatter(s, level)
    case d: DeclarationInterface => formatDeclaration(d, level)
  }

  private def formatCall(call: WdlTaskCall, level: Int): String = {
    val header = s"${highlighter.keyword("call")} ${highlighter.name(call.task.name)}${formatCallAlias(call)}"
    if (call.inputMappings.isEmpty) {
      indent(header, level)
    } else {
      val inputString = call.inputMappings.map {case (k, v) =>
        s"$k=${v.toString(highlighter)}"
      }.mkString(", ")
      indent(s"""$header {
         |  input: $inputString
         |}""".stripMargin, level)
    }
  }

  private def formatScatter(scatter: Scatter, level: Int): String = {
    val children = scatter.children.collect({case c if !c.isInstanceOf[Declaration] => formatScope(c, 1) })
    indent(
      s"""${highlighter.keyword("scatter")} (${scatter.item} in ${scatter.collection.toString(highlighter)}) {
       |${children.mkString("\n")}
       |}""".stripMargin, level)
  }

  private def formatCallAlias(call: WdlTaskCall): String = {
    call.alias.map {a => s" as ${highlighter.alias(a)}"}.getOrElse("")
  }
}
