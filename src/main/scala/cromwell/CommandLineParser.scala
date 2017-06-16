package cromwell

import com.typesafe.config.ConfigFactory
import cromwell.core.path.{DefaultPathBuilder, Path}
import scopt.OptionParser

object CommandLineParser extends App {

  sealed trait Command
  case object Run extends Command
  case object Server extends Command

  case class CommandLineArguments(command: Option[Command] = None,
                                  workflowSource: Option[Path] = None,
                                  workflowInputs: Option[Path] = None,
                                  workflowOptions: Option[Path] = None,
                                  workflowType: Option[String] = Option("WDL"),
                                  workflowTypeVersion: Option[String] = Option("v2.0-draft"),
                                  workflowLabels: Option[Path] = None,
                                  imports: Option[Path] = None,
                                  metadataOutput: Option[Path] = None
                                 )

  lazy val cromwellVersion = ConfigFactory.load("cromwell-version.conf").getConfig("version").getString("cromwell")

  case class ParserAndCommand(parser: OptionParser[CommandLineArguments], command: Option[Command])

  //  cromwell 29
  //  Usage: java -jar /path/to/cromwell.jar [server|run] [options] <args>...
  //
  //    --help                   Cromwell - Workflow Execution Engine
  //    --version
  //  Command: server
  //  Starts a web server on port 8000.  See the web server documentation for more details about the API endpoints.
  //  Command: run [options] workflow-source
  //  Run the workflow and print out the outputs in JSON format.
  //  workflow-source          Workflow source file.
  //  -i, --inputs <value>     Workflow inputs file.
  //  -o, --options <value>    Workflow options file.
  //  -t, --type <value>       Workflow type.
  //  -v, --type-version <value>
  //                           Workflow type version.
  //  -l, --labels <value>     Workflow labels file.
  //  -p, --imports <value>    A directory or zipfile to search for workflow imports.
  //  -m, --metadata-output <value>
  //                           An optional directory path to output metadata.

  def buildParser(): scopt.OptionParser[CommandLineArguments] = {
    new scopt.OptionParser[CommandLineArguments]("java -jar /path/to/cromwell.jar") {
      head("cromwell", cromwellVersion)

      help("help").text("Cromwell - Workflow Execution Engine")

      version("version")

      cmd("server").action((_, c) => c.copy(command = Option(Server))).text(
        "Starts a web server on port 8000.  See the web server documentation for more details about the API endpoints.")

      cmd("run").
        action((_, c) => c.copy(command = Option(Run))).
        text("Run the workflow and print out the outputs in JSON format.").
        children(
          arg[String]("workflow-source").text("Workflow source file.").required().
            action((s, c) => c.copy(workflowSource = Option(DefaultPathBuilder.get(s)))),
          opt[String]('i', "inputs").text("Workflow inputs file.").
            action((s, c) =>
              c.copy(workflowInputs = Option(DefaultPathBuilder.get(s)))),
          opt[String]('o', "options").text("Workflow options file.").
            action((s, c) =>
              c.copy(workflowOptions = Option(DefaultPathBuilder.get(s)))),
          opt[String]('t', "type").text("Workflow type.").
            action((s, c) =>
              c.copy(workflowType = Option(s))),
          opt[String]('v', "type-version").text("Workflow type version.").
            action((s, c) =>
              c.copy(workflowTypeVersion = Option(s))),
          opt[String]('l', "labels").text("Workflow labels file.").
            action((s, c) =>
              c.copy(workflowLabels = Option(DefaultPathBuilder.get(s)))),
          opt[String]('p', "imports").text(
            "A directory or zipfile to search for workflow imports.").
            action((s, c) =>
              c.copy(imports = Option(DefaultPathBuilder.get(s)))),
          opt[String]('m', "metadata-output").text(
            "An optional directory path to output metadata.").
            action((s, c) =>
              c.copy(metadataOutput = Option(DefaultPathBuilder.get(s))))
        )
    }
  }

  def runCromwell(args: CommandLineArguments): Unit = {
    args.command match {
      case Some(Run) => CromwellEntryPoint.runSingle(args)
      case Some(Server) => CromwellEntryPoint.runServer()
      case None => parser.showUsage()
    }
  }

  val parser = buildParser()

  val parsedArgs = parser.parse(args, CommandLineArguments())
  parsedArgs match {
    case Some(pa) => runCromwell(pa)
    case None => parser.showUsage()
  }
}
