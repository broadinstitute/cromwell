package cwl

import ammonite.ops.ImplicitWd._
import ammonite.ops._
import com.typesafe.config.ConfigFactory
import cwl.preprocessor.CwlReference
import org.broadinstitute.heterodon.ExecAndEval

/**
  * Interface for running cwltool.
  */
sealed trait CwltoolRunner {
  def salad(reference: CwlReference): String
}

object CwltoolRunner {
  private lazy val config = ConfigFactory.load

  lazy val instance: CwltoolRunner = {
    val runnerClass = config.getString("cwltool-runner.class")
    Class.forName(runnerClass).getDeclaredConstructor().newInstance().asInstanceOf[CwltoolRunner]
  }
}

/**
  * Runs cwltool as an external process.
  */
final class CwltoolProcess extends CwltoolRunner {
  override def salad(reference: CwlReference): String = {
    val commandResult: CommandResult = %%("cwltool", "--quiet", "--print-pre", reference.pathAsString)
    commandResult.exitCode match {
      case 0 => commandResult.out.string
      case error =>
        throw new RuntimeException(
          s"running CwlTool on file ${reference.pathAsString} resulted in exit code $error and stderr ${commandResult.err.string}")
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

  private def cwltoolSaladEvalStatement(reference: CwlReference): String = s"cwltool_salad('${reference.pathAsString}')"

  def salad(reference: CwlReference): String = {
    val execAndEval = new ExecAndEval()
    execAndEval.apply(cwltoolSaladExecScript, cwltoolSaladEvalStatement(reference)).asInstanceOf[String]
  }
}
