package cromwell

import akka.actor.{Actor, Props}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.services.metadata.{MetadataEvent, MetadataJobKey, MetadataString, MetadataValue}
import cromwell.services.metadata.MetadataService.PutMetadataAction
import MetadataWatchActor._

import scala.concurrent.Promise

// This actor stands in for the service registry and watches for metadata messages that match the optional `matcher`.
// If there is no predicate then use `ignoringBehavior` which ignores all messages.  This is here because there is no
// WorkflowManagerActor in this test that would spin up a real ServiceRegistry.
final case class MetadataWatchActor(promise: Promise[Unit], matchers: Matcher*) extends Actor {

  var unsatisfiedMatchers = matchers

  override def receive = {
    case PutMetadataAction(events) if unsatisfiedMatchers.nonEmpty =>
      unsatisfiedMatchers = unsatisfiedMatchers.filterNot { m => m.matches(events) }
      if (unsatisfiedMatchers.isEmpty) {
        promise.trySuccess(())
        ()
      }
    case PutMetadataAction(_) => // Superfluous message. Ignore
    case _ => throw new Exception("Invalid message to MetadataWatchActor")
  }
}

object MetadataWatchActor {

  def props(promise: Promise[Unit], matchers: Matcher*): Props = Props(MetadataWatchActor(promise, matchers: _*)).withDispatcher(EngineDispatcher)

  trait Matcher {
    private var _fullEventList: List[MetadataEvent] = List.empty
    final def matches(events: Traversable[MetadataEvent]): Boolean = {
      _fullEventList ++= events
      _matches(events)
    }
    def _matches(events: Traversable[MetadataEvent]): Boolean
    private var _nearMisses: List[String] = List.empty
    private def addNearMissInfo(miss: String) = _nearMisses :+= miss
    def nearMissInformation = _nearMisses
    def fullEventList = _fullEventList

    def checkMetadataValueContains(key: String, actual: MetadataValue, expected: String): Boolean = {
      val result = actual.value.contains(expected)
      if (!result) addNearMissInfo(s"Key $key had unexpected value.\nActual value: ${actual.value}\n\nDid not contain: $expected")
      result
    }
  }

  def metadataKeyAttemptChecker(attempt: Int): Option[MetadataJobKey] => Boolean = {
    case Some(jobKey) => jobKey.attempt == attempt
    case None => false
  }

  final case class JobKeyMetadataKeyAndValueContainStringMatcher(jobKeyCheck: Option[MetadataJobKey] => Boolean, key: String, value: String) extends Matcher {
    def _matches(events: Traversable[MetadataEvent]): Boolean = {
      events.exists(e => e.key.key.contains(key) && jobKeyCheck(e.key.jobKey) && e.value.exists { v => v.valueType == MetadataString && checkMetadataValueContains(e.key.key, v, value) })
    }
  }

  abstract class KeyMatchesRegexAndValueContainsStringMatcher(keyTemplate: String, value: String) extends Matcher {
    val templateRegex = keyTemplate.r
    def _matches(events: Traversable[MetadataEvent]): Boolean = {
      events.exists(e => templateRegex.findFirstIn(e.key.key).isDefined &&
        e.value.exists { v => checkMetadataValueContains(e.key.key, v, value) })
    }
  }

  val failurePattern = """failures\[\d*\].message"""
  final case class FailureMatcher(value: String) extends KeyMatchesRegexAndValueContainsStringMatcher(failurePattern, value) { }
}
