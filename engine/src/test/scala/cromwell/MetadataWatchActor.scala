package cromwell

import akka.actor.{Actor, Props}
import cromwell.services.metadata.{MetadataEvent, MetadataJobKey, MetadataString}
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

  def props(promise: Promise[Unit], matchers: Matcher*): Props = Props(MetadataWatchActor(promise, matchers: _*))

  trait Matcher {
    def matches(events: Traversable[MetadataEvent]): Boolean
  }

  def metadataKeyAttemptChecker(attempt: Int): Option[MetadataJobKey] => Boolean = {
    case Some(jobKey) => jobKey.attempt == attempt
    case None => false
  }
  final case class JobKeyMetadataKeyAndValueContainStringMatcher(jobKeyCheck: Option[MetadataJobKey] => Boolean, key: String, value: String) extends Matcher {
    def matches(events: Traversable[MetadataEvent]): Boolean = {
      events.exists(e => e.key.key.contains(key) && jobKeyCheck(e.key.jobKey) && e.value.exists { v => v.valueType == MetadataString && v.value.contains(value) })
    }
  }

  abstract class KeyMatchesRegexAndValueContainsStringMatcher(keyTemplate: String, value: String) extends Matcher {
    val templateRegex = keyTemplate.r
    def matches(events: Traversable[MetadataEvent]): Boolean = {
      events.exists(e => templateRegex.findFirstIn(e.key.key).isDefined && e.value.exists { v => v.value.contains(value) })
    }
  }

  val failurePattern = """failures\[\d*\].message"""
  final case class FailureMatcher(value: String) extends KeyMatchesRegexAndValueContainsStringMatcher(failurePattern, value) {
  }
}
