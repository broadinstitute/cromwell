package womtool

import java.nio.file.Paths

import common.validation.Validation._
import cromwell.core.path.Path
import wdl.draft2.model.formatter.{AnsiSyntaxHighlighter, HtmlSyntaxHighlighter, SyntaxFormatter, SyntaxHighlighter}
import wdl.draft2.model.{AstTools, WdlNamespace}
import womtool.cmdline.HighlightMode.{ConsoleHighlighting, HtmlHighlighting, UnrecognizedHighlightingMode}
import womtool.cmdline._
import womtool.graph.{GraphPrint, WomGraph}
import womtool.input.WomGraphMaker
import womtool.inputs.Inputs
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
    case i: InputsCommandLine => Inputs.inputsJson(i.workflowSource, i.showOptionals)
    case g: WomtoolGraphCommandLine => graph(g.workflowSource.pathAsString)
    case g: WomtoolWomGraphCommandLine => womGraph(g.workflowSource)
    case u: WomtoolWdlV1UpgradeCommandLine => v1upgrade(u.workflowSource.pathAsString)
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

  def parse(workflowSourcePath: String): Termination = {
    SuccessfulTermination(AstTools.getAst(Paths.get(workflowSourcePath)).toPrettyString)
  }

  def v1upgrade(workflowSourcePath: String): Termination = {
    import cats.implicits._
    import common.Checked
    import wdl.draft3.transforms.ast2wdlom.astToFileElement
    import wdl.draft3.transforms.parsing.fileToAst
    import wdl.draft3.transforms.wdlom2wdl.WdlWriter.ops._
    import wdl.draft3.transforms.wdlom2wdl.WdlWriterImpl.fileElementWriter
    import wdl.model.draft3.elements.FileElement

    val loader = fileToAst andThen astToFileElement
    val model: Checked[FileElement] = loader.run(Paths.get(workflowSourcePath))

    model match {
      case Right(wdlModel) => SuccessfulTermination(wdlModel.toWdlV1)
      case Left(errorList) => UnsuccessfulTermination(errorList.toList.mkString("[", ",", "]"))
    }
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

  def womGraph(workflowSourcePath: Path): Termination = {
    WomGraphMaker.fromFiles(mainFile = workflowSourcePath, inputs = None).contextualizeErrors("create wom Graph") match {
      case Right(graph) => SuccessfulTermination (new WomGraph(graphName = "workflow", graph).digraphDot)
      case Left(errors) => UnsuccessfulTermination(errors.toList.mkString(System.lineSeparator, System.lineSeparator, System.lineSeparator))
    }
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
