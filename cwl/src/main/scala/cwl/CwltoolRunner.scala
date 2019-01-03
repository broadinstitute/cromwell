package cwl

import ammonite.ops.ImplicitWd._
import ammonite.ops._
import better.files.{File => BFile}
import com.typesafe.config.ConfigFactory
import org.broadinstitute.heterodon.ExecAndEval

/**
  * Interface for running cwltool.
  */
sealed trait CwltoolRunner {
  def salad(file: BFile): String
}

object CwltoolRunner {
  private lazy val config = ConfigFactory.load

  lazy val instance: CwltoolRunner = {
    val runnerClass = config.getString("cwltool-runner.class")
    Class.forName(runnerClass).newInstance().asInstanceOf[CwltoolRunner]
  }
}

/**
  * Runs cwltool as an external process.
  */
final class CwltoolProcess extends CwltoolRunner {
  override def salad(file: BFile): String = {
    val commandResult: CommandResult = %%("cwltool", "--quiet", "--print-pre", file.toString)
    commandResult.exitCode match {
      case 0 => commandResult.out.string
      case error =>
        throw new RuntimeException(
          s"running CwlTool on file $file resulted in exit code $error and stderr ${commandResult.err.string}")
    }
  }
}

/**
  * Runs cwltool via heterodon.
  *
  * https://github.com/broadinstitute/heterodon
  */
final class CwltoolHeterodon extends CwltoolRunner {
  // Trimmed down version of
  // https://github.com/common-workflow-language/cwltool/blob/1a839255795882894b4bbfea6d909a74cacb1d6a/cwltool/main.py#L355
  private val cwltoolSaladExecScript =
  """|import json
     |import logging
     |
     |from cwltool.load_tool import fetch_document, resolve_tool_uri, validate_document
     |from cwltool.loghandler import _logger
     |
     |
     |def cwltool_salad(path):
     |    _logger.setLevel(logging.WARN)
     |    uri, tool_file_uri = resolve_tool_uri(path)
     |    document_loader, workflowobj, uri = fetch_document(uri)
     |    document_loader, avsc_names, processobj, metadata, uri \
     |        = validate_document(document_loader, workflowobj, uri, preprocess_only=True)
     |    return json.dumps(processobj, indent=4)
     |""".stripMargin

  private def cwltoolSaladEvalStatement(file: BFile): String = s"cwltool_salad('${file.pathAsString}')"

  def salad(file: BFile): String = {
    val execAndEval = new ExecAndEval()
    execAndEval.apply(cwltoolSaladExecScript, cwltoolSaladEvalStatement(file)).asInstanceOf[String]
  }
}
