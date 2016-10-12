package cromwell.backend.impl.htcondor.caching

import akka.actor.{Actor, ActorLogging}
import cromwell.backend.BackendJobExecutionActor.JobSucceededResponse
import cromwell.backend.impl.htcondor.caching.CacheActor._
import cromwell.backend.impl.htcondor.caching.exception.{CachedResultAlreadyExistException, CachedResultNotFoundException}
import cromwell.backend.impl.htcondor.caching.model.CachedExecutionResult

object CacheActor {

  trait CacheActorCommand
  case class ReadExecutionResult(hash: String) extends CacheActorCommand
  case class StoreExecutionResult(hash: String, succeededResponse: JobSucceededResponse) extends CacheActorCommand

  trait CacheActorResponse
  case class ExecutionResultFound(succeededResponse: JobSucceededResponse) extends CacheActorResponse
  case object ExecutionResultNotFound extends CacheActorResponse
  case class ExecutionResultStored(hash: String) extends CacheActorResponse
  case object ExecutionResultAlreadyExist extends CacheActorResponse

}

trait CacheActor extends Actor with ActorLogging {
  def tag: String = "[CacheActor]"
  def forceRewrite: Boolean = false

  override def receive: Receive = {
    case ReadExecutionResult(hash) =>
      try {
        val executionResult = readExecutionResult(hash)
        log.info(s"{} Execution result found in cache for hash {}. Returning result: {}.", tag, hash, executionResult)
        sender() ! ExecutionResultFound(executionResult.succeededResponse)
      } catch {
        case ex: CachedResultNotFoundException => sender() ! ExecutionResultNotFound
      }

    case StoreExecutionResult(hash, succeededResult) =>
      try {
        storeExecutionResult(CachedExecutionResult(hash, succeededResult))
        log.info(s"{} Cache entry for job [{}] stored successfully.", tag, hash)
        sender() ! ExecutionResultStored(hash)
      } catch {
        case ex: CachedResultAlreadyExistException => sender() ! ExecutionResultAlreadyExist
      }
  }

  def readExecutionResult(hash: String): CachedExecutionResult

  def storeExecutionResult(cachedExecutionResult: CachedExecutionResult): Unit

}
