package womtool

import java.nio.file.Paths
import javax.swing.text.Highlighter

import spray.json._
import wdl.draft2.model.formatter.{AnsiSyntaxHighlighter, HtmlSyntaxHighlighter, SyntaxFormatter, SyntaxHighlighter}
import wdl.draft2.model.{AstTools, WdlNamespace, WdlNamespaceWithWorkflow}
import womtool.cmdline.HighlightMode.{ConsoleHighlighting, HtmlHighlighting, UnrecognizedHighlightingMode}
import womtool.cmdline._
import womtool.graph.{GraphPrint, WomGraph}
import womtool.validate.Validate

import scala.util.{Failure, Success}

object WomtoolMain extends App {
  sealed trait Termination {
    def stdout: Option[String]
    def stderr: Option[String]
    def returnCode: Int
  }

  case class SuccessfulTermination(output: String) extends Termination {
    override val stdout: Option[String] = Option(output)
    override val stderr: Option[String] = None
    override val returnCode = 0
  }

  case class UnsuccessfulTermination(output: String) extends Termination {
    override val stdout: Option[String] = None
    override val stderr: Option[String] = Option(output)
    override val returnCode = 1
  }

  case class BadUsageTermination(usage: String) extends Termination {
    override val returnCode = 1
    override val stdout: Option[String] = None
    override val stderr: Option[String] = Option(usage)
  }

  def dispatchCommand(commandLineArgs: ValidatedWomtoolCommandLine): Termination = commandLineArgs match {
    case v: ValidateCommandLine => Validate.validate(v.workflowSource, v.inputs)
    case p: ParseCommandLine => parse(p.workflowSource.pathAsString)
    case h: HighlightCommandLine => highlight(h.workflowSource.pathAsString, h.highlightMode)
    case i: InputsCommandLine => inputs(i.workflowSource.pathAsString)
    case g: WomtoolGraphCommandLine => graph(g.workflowSource.pathAsString)
    case g: WomtoolWomGraphCommandLine => womGraph(g.workflowSource.pathAsString)
    case _ => BadUsageTermination(WomtoolCommandLineParser.instance.usage)
  }


  def highlight(workflowSourcePath: String, mode: HighlightMode): Termination = {

    def highlight(highlighter: SyntaxHighlighter) = {
      loadWdl(workflowSourcePath) { namespace =>
        SuccessfulTermination(new SyntaxFormatter(highlighter).format(namespace))
      }
    }

    mode match {
      case HtmlHighlighting => highlight(HtmlSyntaxHighlighter)
      case ConsoleHighlighting => highlight(AnsiSyntaxHighlighter)
      case UnrecognizedHighlightingMode(m) => BadUsageTermination(s"Unknown highlighter mode: $m")
    }
  }

  def inputs(workflowSourcePath: String): Termination = {
    loadWdl(workflowSourcePath) { namespace =>
      import wom.types.WomTypeJsonFormatter._
      val msg = namespace match {
        case x: WdlNamespaceWithWorkflow => x.workflow.inputs.toJson.prettyPrint
        case _ => "WDL does not have a local workflow"
      }

      SuccessfulTermination(msg)
    }
  }

  def parse(workflowSourcePath: String): Termination = {
    SuccessfulTermination(AstTools.getAst(Paths.get(workflowSourcePath)).toPrettyString)
  }

  def graph(workflowSourcePath: String): Termination = {

      val workflowDigraph = GraphPrint.generateWorkflowDigraph(workflowSourcePath)

      val result = s"""|digraph ${workflowDigraph.workflowName} {
                       |  compound=true;
                       |  ${workflowDigraph.digraph.links.mkString(System.lineSeparator + "  ")}
                       |  ${workflowDigraph.digraph.nodes.mkString(System.lineSeparator + "  ")}
                       |}
                       |"""
      SuccessfulTermination(result.stripMargin)
  }

  def womGraph(workflowSourcePath: String): Termination = {
    SuccessfulTermination(WomGraph.fromFiles(workflowSourcePath).digraphDot)
  }

  private[this] def loadWdl(path: String)(f: WdlNamespace => Termination): Termination = {
    WdlNamespace.loadUsingPath(Paths.get(path), None, None) match {
      case Success(namespace) => f(namespace)
      case Failure(r: RuntimeException) => throw new RuntimeException("Unexpected failure mode", r)
      case Failure(t) => UnsuccessfulTermination(t.getMessage)
    }
  }


  def runWomtool(cmdLineArgs: Seq[String]): Termination = {
    val parsedArgs = WomtoolCommandLineParser.instance.parse(cmdLineArgs, PartialWomtoolCommandLineArguments()) flatMap WomtoolCommandLineParser.validateCommandLine

    parsedArgs match {
      case Some(pa) => dispatchCommand(pa)
      case None => BadUsageTermination(WomtoolCommandLineParser.instance.usage)
    }
  }

  val termination = runWomtool(args)
  termination.stdout foreach Console.out.println
  termination.stderr foreach Console.err.println
  System.exit(termination.returnCode)

}
