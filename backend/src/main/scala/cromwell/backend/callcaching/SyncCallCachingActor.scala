package cromwell.backend.callcaching

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.callcaching.CallCachingActor.{apply => _, _}
import cromwell.backend.callcaching.SyncCallCachingActor._
import wdl4s.RuntimeAttributes

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.language.postfixOps

case class SyncCallCachingActor (
  backendJobDescriptor: BackendJobDescriptor,
  hashFileFunction: String => String,
  hashRuntimeAttributesFunction: RuntimeAttributes => String,
  readFromCache: Boolean,
  writeToCache: Boolean) extends Actor with ActorLogging {

  implicit val ec: ExecutionContext = context.dispatcher

  private var currentClient: Option[ActorRef] = None
  private val callCacheActor = context.actorOf(CallCachingActor.props(backendJobDescriptor, hashFileFunction, hashRuntimeAttributesFunction, readFromCache = readFromCache, writeToCache = writeToCache))
  private val tag = s"SyncCallCachingActor for ${backendJobDescriptor.key.call.unqualifiedName}"

  override def receive: Receive = {
    case CheckCache if currentClient.isEmpty =>
      currentClient = Some(sender)
      callCacheActor ! CheckCache
    case CheckCache =>
      log.error(s"$tag: You shouldn't be using me like this! Only one request at a time please!")
    case CallCacheCheckInProgress =>
      context.system.scheduler.scheduleOnce(500 milliseconds, self, RecheckCacheNow)
    case RecheckCacheNow =>
      callCacheActor ! CheckCache
    case hit: CallCacheHit =>
      currentClient foreach { _ ! hit }
      currentClient = None
    case CallCacheMiss =>
      currentClient foreach { _ ! CallCacheMiss }
      currentClient = None

    case writeCache: WriteCache if currentClient.isEmpty =>
      currentClient = Some(sender)
      callCacheActor ! writeCache
    case WriteCache =>
      log.error(s"$tag: You shouldn't be using me like this! Only one request at a time please!")
    case CallCacheWriteInProgress =>
      context.system.scheduler.scheduleOnce(500 milliseconds, self, RecheckWriteNow)
    case RecheckWriteNow =>
      callCacheActor ! CheckWriteComplete
    case CallCacheWriteSuccess =>
      currentClient foreach { _ ! CallCacheWriteSuccess }
      currentClient = None
    case failure: CallCacheWriteFailure =>
      currentClient foreach { _ ! failure }
      currentClient = None
  }
}

object SyncCallCachingActor {
  def props
  (
    backendJobDescriptor: BackendJobDescriptor,
    hashFileFunction: String => String,
    hashRuntimeAttributesFunction: RuntimeAttributes => String,
    readFromCache: Boolean,
    writeToCache: Boolean) : Props = Props(SyncCallCachingActor(backendJobDescriptor, hashFileFunction, hashRuntimeAttributesFunction, readFromCache = readFromCache, writeToCache = writeToCache))

  private[callcaching] case object RecheckCacheNow
  private[callcaching] case object RecheckWriteNow
}