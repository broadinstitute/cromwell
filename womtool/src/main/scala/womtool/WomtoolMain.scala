package womtool

import java.nio.file.Paths

import com.typesafe.scalalogging.StrictLogging
import common.validation.Validation._
import cromwell.core.path.Path
import cromwell.languages.util.ImportResolver.HttpResolver
import wdl.draft2.model.formatter.{AnsiSyntaxHighlighter, HtmlSyntaxHighlighter, SyntaxFormatter, SyntaxHighlighter}
import wdl.draft2.model.{AstTools, WdlNamespace}
import wom.views.GraphPrint
import womtool.cmdline.HighlightMode.{ConsoleHighlighting, HtmlHighlighting, UnrecognizedHighlightingMode}
import womtool.cmdline._
import womtool.graph.WomGraph
import womtool.input.WomGraphMaker
import womtool.inputs.Inputs
import womtool.outputs.Outputs
import womtool.validate.Validate

import scala.util.{Failure, Success}

object WomtoolMain extends App with StrictLogging {

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
    case v: ValidateCommandLine => Validate.validate(v.workflowSource, v.inputs, v.listDependencies)
    case p: ParseCommandLine => parse(p.workflowSource.pathAsString)
    case h: HighlightCommandLine => highlight(h.workflowSource.pathAsString, h.highlightMode)
    case i: InputsCommandLine => Inputs.inputsJson(i.workflowSource, i.showOptionals)
    case o: OutputsCommandLine => Outputs.outputsJson(o.workflowSource)
    case g: WomtoolGraphCommandLine => graph(g.workflowSource)
    case g: WomtoolWomGraphCommandLine => womGraph(g.workflowSource)
    case _ => BadUsageTermination(WomtoolCommandLineParser.instance.usage)
  }

  def highlight(workflowSourcePath: String, mode: HighlightMode): Termination = {

    def highlight(highlighter: SyntaxHighlighter) =
      loadWdl(workflowSourcePath) { namespace =>
        SuccessfulTermination(new SyntaxFormatter(highlighter).format(namespace))
      }

    mode match {
      case HtmlHighlighting => highlight(HtmlSyntaxHighlighter)
      case ConsoleHighlighting => highlight(AnsiSyntaxHighlighter)
      case UnrecognizedHighlightingMode(m) => BadUsageTermination(s"Unknown highlighter mode: $m")
    }
  }

  def parse(workflowSourcePath: String): Termination =
    SuccessfulTermination(AstTools.getAst(Paths.get(workflowSourcePath)).toPrettyString)

  def graph(workflowSourcePath: Path): Termination =
    WomGraphMaker
      .getBundle(workflowSourcePath)
      .flatMap(_.toExecutableCallable)
      .contextualizeErrors("create wom bundle") match {
      case Right(executable) => SuccessfulTermination(new GraphPrint(executable).dotString)
      case Left(errors) =>
        UnsuccessfulTermination(
          errors.toList.mkString(System.lineSeparator, System.lineSeparator, System.lineSeparator)
        )
    }

  def womGraph(workflowSourcePath: Path): Termination =
    WomGraphMaker
      .fromFiles(mainFile = workflowSourcePath, inputs = None)
      .contextualizeErrors("create wom Graph") match {
      case Right(graphWithImports) =>
        SuccessfulTermination(new WomGraph(graphName = "workflow", graphWithImports.graph).digraphDot)
      case Left(errors) =>
        UnsuccessfulTermination(
          errors.toList.mkString(System.lineSeparator, System.lineSeparator, System.lineSeparator)
        )
    }

  private[this] def loadWdl(path: String)(f: WdlNamespace => Termination): Termination =
    WdlNamespace.loadUsingPath(Paths.get(path), None, None) match {
      case Success(namespace) => f(namespace)
      case Failure(r: RuntimeException) => throw new RuntimeException("Unexpected failure mode", r)
      case Failure(t) => UnsuccessfulTermination(t.getMessage)
    }

  def runWomtool(cmdLineArgs: Seq[String]): Termination = {
    val parsedArgs = WomtoolCommandLineParser.instance.parse(cmdLineArgs,
                                                             PartialWomtoolCommandLineArguments()
    ) flatMap WomtoolCommandLineParser.validateCommandLine

    parsedArgs match {
      case Some(pa) => dispatchCommand(pa)
      case None => BadUsageTermination(WomtoolCommandLineParser.instance.usage)
    }
  }

  val termination = runWomtool(args.toIndexedSeq)
  termination.stdout foreach Console.out.println
  termination.stderr foreach Console.err.println

  HttpResolver.closeBackendIfNecessary()

  System.exit(termination.returnCode)

}
