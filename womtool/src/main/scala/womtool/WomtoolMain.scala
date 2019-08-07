package womtool

import java.nio.file.Paths

import better.files.File
import common.validation.Validation._
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.languages.util.ImportResolver.HttpResolver
import languages.wdl.draft2.WdlDraft2LanguageFactory
import wdl.draft2.model.formatter.{AnsiSyntaxHighlighter, HtmlSyntaxHighlighter, SyntaxFormatter, SyntaxHighlighter}
import wdl.draft2.model.{AstTools, WdlNamespace}
import wdl.transforms.base.wdlom2wdl.WdlWriter.ops._
import wdl.transforms.base.wdlom2wdl.WdlWriterImpl.fileElementWriter
import womtool.cmdline.HighlightMode.{ConsoleHighlighting, HtmlHighlighting, UnrecognizedHighlightingMode}
import womtool.cmdline._
import womtool.graph.{GraphPrint, WomGraph}
import womtool.input.WomGraphMaker
import womtool.inputs.Inputs
import womtool.validate.Validate
import womtool.wom2wdlom.WomToWdlom.womBundleToFileElement

import scala.util.{Failure, Success, Try}

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
    case v: ValidateCommandLine => Validate.validate(v.workflowSource, v.inputs, v.listDependencies)
    case p: ParseCommandLine => parse(p.workflowSource.pathAsString)
    case h: HighlightCommandLine => highlight(h.workflowSource.pathAsString, h.highlightMode)
    case i: InputsCommandLine => Inputs.inputsJson(i.workflowSource, i.showOptionals)
    case g: WomtoolGraphCommandLine => graph(g.workflowSource.pathAsString)
    case g: WomtoolWomGraphCommandLine => womGraph(g.workflowSource)
    case u: WomtoolWdlUpgradeCommandLine => upgrade(u.workflowSource.pathAsString)
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

  def upgrade(workflowSourcePath: String): Termination = {
    import wdl.draft2.model.Import
    import wdl.model.draft3.elements.ImportElement

    // Get imports directly from WdlNamespace, because they are erased during WOMification
    val maybeWdlNamespace: Try[WdlNamespace] =
      WdlNamespace.loadUsingPath(
        Paths.get(workflowSourcePath),
        None,
        Option(List(
          WdlNamespace.directoryResolver(File(workflowSourcePath).parent),
          WdlNamespace.fileResolver,
          WdlDraft2LanguageFactory.httpResolver
        ))
      )

    def upgradeImport(draft2Import: Import): ImportElement = {
      if (draft2Import.namespaceName.nonEmpty)
        // draft-2 does not have structs, so the source WDL will not have any for us to rename
        ImportElement(draft2Import.uri, Option(draft2Import.namespaceName), Map())
      else
        ImportElement(draft2Import.uri, None, Map())
    }

    val maybeWdl: Try[Path] = DefaultPathBuilder.build(workflowSourcePath)

    (maybeWdl, maybeWdlNamespace) match {
      case (Success(wdl), Success(wdlNamespace)) =>
        val maybeWomBundle = WomGraphMaker.getBundle(wdl)
        maybeWomBundle match {
          case Right(womBundle) =>
            val maybeFileElement = womBundleToFileElement.run(womBundle)
            maybeFileElement match {
              case Right(fileElement) =>
                SuccessfulTermination(
                  fileElement.copy(imports = wdlNamespace.imports.map(upgradeImport)).toWdlV1)
              case Left(errors) =>
                UnsuccessfulTermination(s"WDL parsing succeeded but could not create WOM: ${errors.toList.mkString("[", ",", "]")}")
            }
          case Left(errors) =>
            UnsuccessfulTermination(s"WDL parsing succeeded but could not create WOM: ${errors.toList.mkString("[", ",", "]")}")
        }
      case (Failure(throwable), _) =>
        UnsuccessfulTermination(s"Failed to load WDL source: ${throwable.getMessage}")
      case (_, Failure(throwable)) =>
        UnsuccessfulTermination(s"Failed to load WDL source: ${throwable.getMessage}")
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
      case Right(graphWithImports) => SuccessfulTermination (new WomGraph(graphName = "workflow", graphWithImports.graph).digraphDot)
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

  HttpResolver.closeBackendIfNecessary()

  System.exit(termination.returnCode)

}
