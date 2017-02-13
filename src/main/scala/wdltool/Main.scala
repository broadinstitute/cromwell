package wdltool

import java.nio.file.Paths

import wdl4s.formatter.{AnsiSyntaxHighlighter, HtmlSyntaxHighlighter, SyntaxFormatter}
import wdl4s._
import spray.json._


import scala.util.{Failure, Success, Try}

object Main extends App {
  sealed trait Termination {
    def output: String
    def returnCode: Int
  }

  case class SuccessfulTermination(output: String) extends Termination {
    override val returnCode = 0
  }

  case class UnsuccessfulTermination(output: String) extends Termination {
    override val returnCode = 1
  }

  case object BadUsageTermination extends Termination {
    override val returnCode = 1
    override val output = UsageMessage
  }

  def dispatchCommand(args: Seq[String]): Termination = {
    getAction(args) match {
      case Some(x) if x == Actions.Validate => validate(args.tail)
      case Some(x) if x == Actions.Highlight => highlight(args.tail)
      case Some(x) if x == Actions.Inputs => inputs(args.tail)
      case Some(x) if x == Actions.Parse => parse(args.tail)
      case Some(x) if x == Actions.Graph => graph(args.tail)
      case _ => BadUsageTermination
    }
  }

  def validate(args: Seq[String]): Termination = {
    continueIf(args.length == 1) {
      loadWdl(args.head) { _ => SuccessfulTermination("") }
    }
  }

  def highlight(args: Seq[String]): Termination = {
    continueIf(args.length == 2 && Seq("html", "console").contains(args(1))) {
      loadWdl(args.head) { namespace =>
        val formatter = new SyntaxFormatter(if (args(1) == "html") HtmlSyntaxHighlighter else AnsiSyntaxHighlighter)
        SuccessfulTermination(formatter.format(namespace))
      }
    }
  }

  def inputs(args: Seq[String]): Termination = {
    continueIf(args.length == 1) {
      loadWdl(args.head) { namespace =>
        import wdl4s.types.WdlTypeJsonFormatter._
        val msg = namespace match {
          case x: WdlNamespaceWithWorkflow => x.workflow.inputs.toJson.prettyPrint
          case _ => "WDL does not have a local workflow"
        }

        SuccessfulTermination(msg)
      }
    }
  }

  def parse(args: Seq[String]): Termination = {
    continueIf(args.length == 1) {
      SuccessfulTermination(AstTools.getAst(Paths.get(args.head)).toPrettyString)
    }
  }

  def graph(args: Seq[String]): Termination = {
    continueIf(args.length == 1 || (args.length == 2 && args.head.equals("--all"))) {

      val (file, allNodesMode) =
        if (args.size == 1) (args.head, false)
        else (args(1), true)

      val workflowDigraph = GraphPrint.generateWorkflowDigraph(file, allNodesMode)

      val result = s"""|digraph ${workflowDigraph.workflowName} {
                       |  ${workflowDigraph.digraph.mkString(System.lineSeparator + "  ")}
                       |}
                       |"""
      SuccessfulTermination(result.stripMargin)
    }
  }

  private[this] def continueIf(valid: => Boolean)(block: => Termination): Termination = if (valid) block else BadUsageTermination

  private[this] def loadWdl(path: String)(f: WdlNamespace => Termination): Termination = {
    WdlNamespace.loadUsingPath(Paths.get(path), None, None) match {
      case Success(namespace) => f(namespace)
      case Failure(t) => UnsuccessfulTermination(t.getMessage)
    }
  }

  private def getAction(args: Seq[String]): Option[Actions.Value] = for {
    arg <- args.headOption
    argCapitalized = arg.capitalize
    action <- Actions.values find (_.toString == argCapitalized)
  } yield action

  object Actions extends Enumeration {
    val Parse, Validate, Highlight, Inputs, Graph = Value
  }

  val UsageMessage = """
                       |java -jar wdltool.jar <action> <parameters>
                       |
                       |Actions:
                       |validate <WDL file>
                       |
                       |  Performs full validation of the WDL file including syntax
                       |  and semantic checking
                       |
                       |inputs <WDL file>
                       |
                       |  Print a JSON skeleton file of the inputs needed for this
                       |  workflow.  Fill in the values in this JSON document and
                       |  pass it in to the 'run' subcommand.
                       |
                       |highlight <WDL file> <html|console>
                       |
                       |  Reformats and colorizes/tags a WDL file. The second
                       |  parameter is the output type.  "html" will output the WDL
                       |  file with <span> tags around elements.  "console" mode
                       |  will output colorized text to the terminal
                       |
                       |parse <WDL file>
                       |
                       |  Compares a WDL file against the grammar and prints out an
                       |  abstract syntax tree if it is valid, and a syntax error
                       |  otherwise.  Note that higher-level AST checks are not done
                       |  via this sub-command and the 'validate' subcommand should
                       |  be used for full validation.
                       |graph [--all] <WDL file>
                       |
                       |  Reads a WDL file against the grammar and prints out a
                       |  .dot of the DAG if it is valid, and a syntax error
                       |  otherwise.
                       |  Use [--all] to show all graph nodes in the WDL spec,
                       |  even the non-executable nodes.
                     """.stripMargin

  val termination = dispatchCommand(args)

  termination match {
    case SuccessfulTermination(s) => println(s)
    case UnsuccessfulTermination(s) => Console.err.println(s)
    case BadUsageTermination => Console.err.println(UsageMessage)
  }

  termination.returnCode
}
