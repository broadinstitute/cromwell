package centaur.test.metadata

import java.time.OffsetDateTime

import cats.effect._
import cats.instances.either._
import cats.instances.option._
import cats.instances.vector._
import cats.syntax.apply._
import cats.syntax.traverse._
import io.circe._
import io.circe.parser._

/**
  * The first failure message of a call. Based on the JMUI FailureMessage:
  *
  * https://github.com/DataBiosphere/job-manager/blob/f83e4284e2419389b7e515720c9d960d2eb81a29/servers/cromwell/jobs/controllers/jobs_controller.py#L155-L162
  */
case class CallAttemptFailure
(
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
  def buildFailures(jsonOption: Option[String]): IO[Vector[CallAttemptFailure]] = {
    jsonOption.map(buildFailures).getOrElse(IO.pure(Vector.empty))
  }

  def buildFailures(json: String): IO[Vector[CallAttemptFailure]] = {
    IO.fromEither(decode[Vector[CallAttemptFailure]](json))
  }

  private implicit val decodeFailures: Decoder[Vector[CallAttemptFailure]] = {
    Decoder.instance { c =>
      for {
        workflowId <- c.get[String]("id")
        calls <- c.get[Map[String, Json]]("calls").map(_.toVector)
        callAttemptFailures <- calls.flatTraverse[Decoder.Result, CallAttemptFailure] {
          case (callName, callJson) =>
            val decoderCallAttempt = decodeFromCallAttempt(workflowId, callName)
            callJson.as[Vector[Option[CallAttemptFailure]]](Decoder.decodeVector(decoderCallAttempt)).map(_.flatten)
        }
      } yield callAttemptFailures
    } or Decoder.const(Vector.empty)
  }

  private def decodeFromCallAttempt(workflowId: String, callName: String): Decoder[Option[CallAttemptFailure]] = {
    Decoder.instance { c =>
      for {
        shardIndexOption <- c.get[Option[Int]]("shardIndex")
        attemptOption <- c.get[Option[Int]]("attempt")
        messageOption <- c.downField("failures").downArray.get[Option[String]]("message")
        startOption <- c.get[Option[OffsetDateTime]]("start")
        endOption <- c.get[Option[OffsetDateTime]]("end")
        stdoutOption <- c.get[Option[String]]("stdout")
        stderrOption <- c.get[Option[String]]("stderr")
        callRootOption <- c.get[Option[String]]("callRoot")
        callAttemptFailureOption = (shardIndexOption, attemptOption, messageOption) mapN {
          (shardIndex, attempt, message) =>
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
      } yield callAttemptFailureOption
    } or Decoder.const(None)
  }
}
