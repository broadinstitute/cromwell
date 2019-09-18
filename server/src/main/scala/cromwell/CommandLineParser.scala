package cromwell

import java.net.URL

import common.util.VersionUtil
import cromwell.CommandLineParser._
import cromwell.CromwellApp._
import cromwell.core.path.DefaultPathBuilder

object CommandLineParser {
  lazy val cromwellVersion = VersionUtil.getVersion("cromwell")
}

//  cromwell 29
//  Usage: java -jar /path/to/cromwell.jar [server|run] [options] <args>...
//
//    --help                   Cromwell - Workflow Execution Engine
//    --version
//  Command: server
//  Starts a web server on port 8000.  See the web server documentation for more details about the API endpoints.
//  Command: run [options] workflow-source
//  Run the workflow and print out the outputs in JSON format.
//  workflow-source          Workflow source file or workflow url.
//  --workflow-root <value>  Workflow root
//  -i, --inputs <value>     Workflow inputs file.
//  -o, --options <value>    Workflow options file.
//  -t, --type <value>       Workflow type.
//  -v, --type-version <value>
//                           Workflow type version.
//  -l, --labels <value>     Workflow labels file.
//  -p, --imports <value>    A directory or zipfile to search for workflow imports.
//  -m, --metadata-output <value>
//                           An optional JSON file path to output metadata.
//  -h, --host               Cromwell server URL
class CommandLineParser extends scopt.OptionParser[CommandLineArguments]("java -jar /path/to/cromwell.jar") {
  
  private def commonSubmissionArguments = List(
    arg[String]("workflow-source").text("Workflow source file or workflow url.").required().
      action((s, c) =>
        c.copy(workflowSource = Option(s))),
    opt[String]("workflow-root").text("Workflow root.").
      action((s, c) =>
        c.copy(workflowRoot = Option(s))),
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
        c.copy(imports = Option(DefaultPathBuilder.get(s))))
  )

  head("cromwell", cromwellVersion)

  help("help").text("Cromwell - Workflow Execution Engine")

  version("version")

  cmd("server").action((_, c) => c.copy(command = Option(Server))).text(
    "Starts a web server on port 8000.  See the web server documentation for more details about the API endpoints.")

  cmd("run").
    action((_, c) => c.copy(command = Option(Run))).
    text("Run the workflow and print out the outputs in JSON format.").
    children(
      commonSubmissionArguments ++ List(
        opt[String]('m', "metadata-output").text(
          "An optional JSON file path to output metadata.").
          action((s, c) =>
            c.copy(metadataOutput = Option(DefaultPathBuilder.get(s))))
      ): _*
    )

  cmd("submit")
    .action((_, c) => c.copy(command = Option(Submit))).
    text("Submit the workflow to a Cromwell server.").
    children(
      commonSubmissionArguments ++ List(
        opt[String]('h', "host").text("Cromwell server URL.").
          action((h, c) =>
            c.copy(host = new URL(h)))
      ): _*
    )
}
