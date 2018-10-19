package cromwell.webservice

import cats.effect.IO
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import cromwell.webservice.metadata.StreamMetadataBuilder
import org.reactivestreams.{Publisher, Subscriber, Subscription}
import org.scalatest.Assertion
import spray.json._

import scala.concurrent.Future

/**
  * To re-use the existing Metadata tests in MetadataBuilderActorSpec, simply override the assertion methods to use the
  * streaming mechanism, keeping all the same test cases
  */
class StreamMetadataBuilderSpec extends MetadataBuilderActorSpec("StreamMetadataBuilder") {
  
  override def assertMetadataResponse(action: MetadataServiceAction,
                             queryReply: MetadataQuery,
                             events: Seq[MetadataEvent],
                             expectedRes: String): Future[Assertion] = {
    val publisher = new TestPublisher(events)
    val metadataStreamer = new StreamMetadataBuilder(defaultTimeout) {
      override def queryToPublisher(query: MetadataQuery) = IO.pure(publisher)
    }
    val result = action match {
      case singleAction: GetSingleWorkflowMetadataAction => metadataStreamer.workflowMetadataQuery(singleAction)
      case queryAction: GetMetadataQueryAction => metadataStreamer.workflowMetadataQuery(queryAction.key)
      case _ => fail(s"Unstreamable action: ${action.getClass.getSimpleName}")
    }
    result.unsafeToFuture() map { b => b shouldBe expectedRes.parseJson}
  }
}

// A publisher that distributes a precomputed sequence of events
class TestPublisher(events: Seq[MetadataEvent]) extends Publisher[MetadataEvent] {
  override def subscribe(s: Subscriber[_ >: MetadataEvent]) = {
    s.onSubscribe( new Subscription {
      private var index: Int = 0
      override def request(n: Long) = {
        (0 until n.toInt).foreach(_ => next())
      }
      
      override def cancel() = s.onComplete()
      
      // Not the most efficient but good enough for tests
      private def next() = {
        if (index < events.length) s.onNext(events(index))
        if (index == events.length - 1) s.onComplete()
        index += 1
      }
    })
  }
}
