package centaur.test.metadata

import java.time.OffsetDateTime

import cats.effect.IO
import cats.implicits._
import cromwell.api.model.OffsetDateTimeJsonFormat
import spray.json._
import spray.json.DefaultJsonProtocol._

import scala.util.Try

/**
  * The first failure message of a call. Based on the JMUI FailureMessage:
  *
  * https://github.com/DataBiosphere/job-manager/blob/f83e4284e2419389b7e515720c9d960d2eb81a29/servers/cromwell/jobs/controllers/jobs_controller.py#L155-L162
  */
case class CallAttemptFailure(
  workflowId: String,
  callFullyQualifiedName: String,
  jobIndex: Int,
  jobAttempt: Int,
  message: String,
  startOption: Option[OffsetDateTime],
  endOption: Option[OffsetDateTime],
  stdoutOption: Option[String],
  stderrOption: Option[String],
  callRootOption: Option[String]
)

object CallAttemptFailure {
  def buildFailures(jsonOption: Option[String]): IO[Vector[CallAttemptFailure]] =
    jsonOption.map(buildFailures).getOrElse(IO.pure(Vector.empty))

  def buildFailures(json: String): IO[Vector[CallAttemptFailure]] =
    IO.fromTry(Try(json.parseJson)) map { parsed =>
      Try(decodeFailures(parsed)).getOrElse(Vector.empty)
    }

  private def decodeFailures(json: JsValue): Vector[CallAttemptFailure] = {
    val fields = json.asJsObject.fields
    val workflowId = fields("id").convertTo[String]
    val calls = fields("calls").asJsObject.fields.toVector
    calls flatMap { case (callName, callJson) =>
      callJson.convertTo[Vector[JsValue]] flatMap { attempt =>
        Try(decodeFromCallAttempt(workflowId, callName, attempt.asJsObject)).toOption.flatten
      }
    }
  }

  private def field[A: JsonReader](json: JsObject, name: String): Option[A] =
    json.fields.get(name).filter(_ != JsNull).map(_.convertTo[A])

  private def decodeFromCallAttempt(workflowId: String, callName: String, json: JsObject): Option[CallAttemptFailure] = {
    val shardIndexOption = field[Int](json, "shardIndex")
    val attemptOption = field[Int](json, "attempt")
    val messageOption = json.fields.get("failures").collect { case JsArray(failures) if failures.nonEmpty =>
      field[String](failures.head.asJsObject, "message")
    }.flatten
    val startOption = field[OffsetDateTime](json, "start")
    val endOption = field[OffsetDateTime](json, "end")
    val stdoutOption = field[String](json, "stdout")
    val stderrOption = field[String](json, "stderr")
    val callRootOption = field[String](json, "callRoot")
    (shardIndexOption, attemptOption, messageOption) mapN { (shardIndex, attempt, message) =>
      new CallAttemptFailure(
        workflowId = workflowId,
        callFullyQualifiedName = callName,
        jobIndex = shardIndex,
        jobAttempt = attempt,
        message = message,
        startOption = startOption,
        endOption = endOption,
        stdoutOption = stdoutOption,
        stderrOption = stderrOption,
        callRootOption = callRootOption
      )
    }
  }
}
