package cwl

import java.util.concurrent.Executors

import ammonite.ops.ImplicitWd._
import ammonite.ops._
import cats.effect.IO
import com.typesafe.config.{Config, ConfigFactory}
import cwl.preprocessor.{CwlFileReference, CwlHttpReference, CwlReference}
import io.circe.Json
import net.ceedubs.ficus.Ficus._
import org.broadinstitute.heterodon.ExecAndEval
import org.http4s.circe._
import org.http4s.client.blaze.BlazeClientBuilder
import org.http4s.{EntityEncoder, Method, ParseFailure, Request, Uri}

import scala.concurrent.ExecutionContext

/**
  * Interface for running cwltool.
  */
sealed trait CwltoolRunner {
  def salad(reference: CwlReference): IO[String]
}

object CwltoolRunner {
  private lazy val config = ConfigFactory.load

  lazy val instance: CwltoolRunner = {
    val runnerClass = config.getString("cwltool-runner.class")
    Class.forName(runnerClass).newInstance().asInstanceOf[CwltoolRunner]
  }
  
  lazy val instanceConfig = config.getAs[Config]("cwltool-runner.config")
}

/**
  * Runs cwltool as an external process.
  */
final class CwltoolProcess extends CwltoolRunner {
  override def salad(reference: CwlReference): IO[String] = {
    val commandResult: IO[CommandResult] = IO(%%("cwltool", "--quiet", "--print-pre", reference.pathAsString))

    commandResult.flatMap({ result =>
      result.exitCode match {
        case 0 => IO.pure(result.out.string)
        case error =>
          IO.raiseError(new RuntimeException(
            s"running CwlTool on file ${reference.pathAsString} resulted in exit code $error and stderr ${result.err.string}"))
      }
    })
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

  def salad(reference: CwlReference): IO[String] = {
    val execAndEval = new ExecAndEval()
    IO(execAndEval.apply(cwltoolSaladExecScript, cwltoolSaladEvalStatement(reference)).asInstanceOf[String])
  }
}

final class CwltoolHttp extends CwltoolRunner {
  private val ec = ExecutionContext.fromExecutor(Executors.newCachedThreadPool())
  private implicit val cs = IO.contextShift(ec)
  private val (client, _) = BlazeClientBuilder[IO](ec).allocate.unsafeRunSync()
  // Uri to send the requests to
  private val uri: Uri = CwltoolRunner.instanceConfig
    .map(_.as[String]("url"))
    .map(Uri.fromString) match {
    case Some(Right(url)) => url
    case Some(Left(f: ParseFailure)) => throw new RuntimeException(s"Invalid url for CwltoolHttp: ${f.message}")
    case None => throw new RuntimeException("Cannot use CwltoolHttp without an URL configured at cwltool-runner.config.url")
  }
  
  override def salad(reference: CwlReference): IO[String] = {
    val request = reference match {
      case CwlFileReference(file, _) => makeRequest(file.contentAsString)
      case CwlHttpReference(pathAsString, _) => makeRequest(Json.obj("path" -> Json.fromString(pathAsString)))
    }

    client.expect[String](request)
  }
  
  private def makeRequest[A](entity: A)(implicit encoder: EntityEncoder[IO, A]) = {
    Request[IO](
      method = Method.POST,
      uri = uri
    ).withEntity(entity)
  }
}
