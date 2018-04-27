package womtool

import java.nio.file.Paths

import cats.data.NonEmptyList
import common.Checked
import common.transforms.CheckedAtoB
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
    case u: WomtoolDraft3UpgradeCommandLine => d3upgrade(u.workflowSource.pathAsString)
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

  // TODO: move to an appropriate location
  import wdl.draft2.parser.WdlParser.Ast
  import wdl.model.draft3.elements.FileElement
  import wdl.model.draft3.elements.TaskDefinitionElement
//  import wdl.model.draft3.elements.CommandSectionElement
//  import wdl.model.draft3.elements.FileBodyElement
  import wdl.model.draft3.elements.ImportElement
  import wdl.draft2.model.AstTools._
  import wdl.model.draft3.elements.WorkflowDefinitionElement
  import wdl.model.draft3.elements.InputsSectionElement
  import wdl.model.draft3.elements.OutputsSectionElement
  import wdl.model.draft3.elements.{OutputDeclarationElement, PrimitiveTypeElement}
  import wdl.model.draft3.elements.TypeElement
//  import wdl.draft2.parser.WdlParser.AstNode
  import wom.types.WomStringType
  import wdl.model.draft3.elements.ExpressionElement.StringLiteral

  def convertAstToFile(ast: Ast): Checked[FileElement] = ast.getName match {
    case "Namespace" =>
      val bodyElements: Seq[Ast] = ast.getAttribute("body").astListAsVector.collect({ case a: Ast => a })
      val imports: Seq[Ast] = ast.getAttribute("imports").astListAsVector.collect({ case a: Ast => a })

      val tasks = bodyElements.filter(_.getName == "Task")
      val workflows = bodyElements.filter(_.getName == "Workflow")

      // TODO: better way to provide this guarantee - match statement with _ that throws exception?
      assert(tasks.size + workflows.size == bodyElements.size, "Missed something!")

      Right(FileElement(
        workflows = workflows.map(convertAstToWorkflow),
        imports = imports.map(convertAstToImport),
        structs = Seq(), // no structs in draft-2 I think?
        tasks = tasks.map(convertAstToTask)))
    case _ => Left(NonEmptyList("What are you doing?", Nil))
  }

  def convertAstToWorkflow(ast: Ast): WorkflowDefinitionElement = {
    val name = ast.getAttribute("name").sourceString
    val bodyElements: Seq[Ast] = ast.getAttribute("body").astListAsVector.collect({ case a: Ast => a })

    val outputs: Ast = bodyElements.filter(_.getName == "WorkflowOutputs").head

    WorkflowDefinitionElement(
      name = name,
      inputsSection = None,
      graphElements = Set.empty,
      outputsSection = convertAstToOutputsSection(outputs),
      metaSection = None,
      parameterMetaSection = None
    )
  }

  def convertAstToOutputsSection(ast: Ast): Option[OutputsSectionElement] = {
    val outputElements: Seq[Ast] = ast.getAttribute("outputs").astListAsVector.collect({ case a: Ast => a })

    if (outputElements.nonEmpty)
      Some(OutputsSectionElement(outputs = outputElements.map(convertAstToOutputDeclarationElement)))
    else
      None
  }

  def convertAstToOutputDeclarationElement(ast: Ast): OutputDeclarationElement = {
    //    val typeElement = ast.getAttribute("typeElement")
    //    val name = ast.getAttribute("name")
    //    val expression = ast.getAttribute("expression")

    OutputDeclarationElement(
      typeElement = PrimitiveTypeElement(WomStringType),
      name = "asdf",
      expression = StringLiteral("wasd")
    )
  }

  def convertAstToTypeElement(ast: Ast): TypeElement = ???

  def convertAstToInputsSectionElement(ast: Ast): InputsSectionElement = ???

  def convertAstToTask(ast: Ast): TaskDefinitionElement = ???

  def convertAstToImport(ast: Ast): ImportElement = ???

  def convertFileElementToString(fileModel: FileElement): Checked[String] = Right("Meaningless string!")

//  object AstToFileBodyElement {
//    def convert(ast: Ast): Checked[FileBodyElement] = ast.getName match {
//      case "Workflow" => astNodeToWorkflowDefinitionElement(ast)
//      case "Task" => astNodeToTaskDefinitionElement(ast)
//      case "Struct" => astNodeToStructEntry(ast)
//      case other => s"No conversion defined for Ast with name $other".invalidNelCheck
//    }
//  }

  def d3upgrade(workflowSourcePath: String): Termination = {
    import cats.implicits._

    def astToModelConverter: CheckedAtoB[Ast, FileElement] = CheckedAtoB.fromCheck(convertAstToFile)
    def modelToStringConverter: CheckedAtoB[FileElement, String] = CheckedAtoB.fromCheck(convertFileElementToString)

    val ast: wdl.draft2.parser.WdlParser.Ast = AstTools.getAst(Paths.get(workflowSourcePath))
    val converted: Checked[String] = (astToModelConverter andThen modelToStringConverter).run(ast)

    converted match {
      case Right(wdl) => SuccessfulTermination(wdl)
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
