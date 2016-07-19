package cromwell.backend

import java.util.UUID

import scala.concurrent.{ExecutionContext, Future}

object ExecutionHash {
  // TODO: PBE: ideally hashes should be deterministic
  def completelyRandomExecutionHash(implicit ec: ExecutionContext): Future[ExecutionHash] = Future.successful(
    ExecutionHash(UUID.randomUUID().toString, dockerHash = None))
}

final case class ExecutionHash(overallHash: String, dockerHash: Option[String])
