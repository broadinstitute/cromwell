package cromwell.services.metadata.impl.pubsub

import java.time.OffsetDateTime

import akka.actor.{ActorInitializationException, Props}
import akka.testkit.EventFilter
import com.google.api.client.auth.oauth2.Credential
import com.google.api.services.pubsub.model.Topic
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.WorkflowId
import cromwell.services.ServicesSpec
import cromwell.services.metadata.MetadataService.{MetadataWriteFailure, MetadataWriteSuccess, PutMetadataAction, PutMetadataActionAndRespond}
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import org.broadinstitute.dsde.workbench.google.GooglePubSubDAO
import org.broadinstitute.dsde.workbench.google.GooglePubSubDAO.PubSubMessage

import scala.concurrent.{ExecutionContext, Future}

class PubSubMetadataServiceActorSpec extends ServicesSpec("PubSubMetadata") {
  import PubSubMetadataServiceActorSpec._

  "A PubSubMetadataActor with an empty serviceConfig" should {
    "fail to build" in {
      EventFilter[ActorInitializationException](occurrences = 1) intercept {
        system.actorOf(Props(new PubSubMetadataServiceActor(emptyConfig, emptyConfig)))
      }
    }
  }

  "A PubSubMetadataActor with a subscription" should {
    "should ensure topic exists" in {
      EventFilter.info("Ensuring topic bar exists", occurrences = 1) intercept {
        system.actorOf(Props(new SuccessfulMockPubSubMetadataServiceActor(configWithSubscription, emptyConfig)))
      }
    }

    "should create the requested subscription" in {
      EventFilter.info("Creating subscription baz", occurrences = 1) intercept {
        system.actorOf(Props(new SuccessfulMockPubSubMetadataServiceActor(configWithSubscription, emptyConfig)))
      }
    }

    "should log to debug on a publish request" in {
      EventFilter.debug(start = "Publishing to", occurrences = 1) intercept {
        val psma = system.actorOf(Props(new SuccessfulMockPubSubMetadataServiceActor(configWithSubscription, emptyConfig)))
        psma ! PutMetadataAction(List(Event))
      }
    }

    "should ack on a publish-ack request" in {
      val psma = system.actorOf(Props(new SuccessfulMockPubSubMetadataServiceActor(configWithSubscription, emptyConfig)))
      psma ! PutMetadataActionAndRespond(List(Event), testActor)
      expectMsgClass(classOf[MetadataWriteSuccess])
    }
  }

  "A PubSubMetadataActor without a subscription" should {
    "should ensure topic exists" in {
      EventFilter.info("Ensuring topic bar exists", occurrences = 1) intercept {
        system.actorOf(Props(new SuccessfulMockPubSubMetadataServiceActor(configWithoutSubscription, emptyConfig)))
      }
    }

    "should not be creating a subscription" in {
      EventFilter.info("Not creating a subscription", occurrences = 1) intercept {
        system.actorOf(Props(new SuccessfulMockPubSubMetadataServiceActor(configWithoutSubscription, emptyConfig)))
      }
    }

    "should log to debug on a publish request" in {
      EventFilter.debug(start = "Publishing to", occurrences = 1) intercept {
        val psma = system.actorOf(Props(new SuccessfulMockPubSubMetadataServiceActor(configWithoutSubscription, emptyConfig)))
        psma ! PutMetadataAction(List(Event))
      }
    }

    "should ack on a publish-ack request" in {
      val psma = system.actorOf(Props(new SuccessfulMockPubSubMetadataServiceActor(configWithoutSubscription, emptyConfig)))
      psma ! PutMetadataActionAndRespond(List(Event), testActor)
      expectMsgClass(classOf[MetadataWriteSuccess])
    }
  }

  "A PubSubMetadataActor who fails at all things" should {
    "fail to create a topic" in {
      EventFilter[RuntimeException](start = "Unable to create topic", occurrences = 1) intercept {
        system.actorOf(Props(new FailingToCreateTopicMockPubSubMetadataServiceActor(configWithSubscription, emptyConfig)))
      }
    }
  }

  "A PubSubMetadataActor who fails to publish" should {
    "should ensure topic exists" in {
      EventFilter.info("Ensuring topic bar exists", occurrences = 1) intercept {
        system.actorOf(Props(new FailToPublishMockPubSubMetadataServiceActor(configWithoutSubscription, emptyConfig)))
      }
    }

    "should log to error on a publish request" in {
      EventFilter[RuntimeException](start = "Failed to post metadata: ", occurrences = 1) intercept {
        val psma = system.actorOf(Props(new FailToPublishMockPubSubMetadataServiceActor(configWithoutSubscription, emptyConfig)))
        psma ! PutMetadataAction(List(Event))
      }
    }

    "should ack on a publish-ack request" in {
      val psma = system.actorOf(Props(new FailToPublishMockPubSubMetadataServiceActor(configWithoutSubscription, emptyConfig)))
      psma ! PutMetadataActionAndRespond(List(Event), testActor)
      expectMsgClass(classOf[MetadataWriteFailure])
    }
  }
}

object PubSubMetadataServiceActorSpec {
  /** A variant of PubSubMetadataServiceActor with a GooglePubSubDAO which will always return success */
  class SuccessfulMockPubSubMetadataServiceActor(serviceConfig: Config, globalConfig: Config)
    extends PubSubMetadataServiceActor(serviceConfig, globalConfig) {

    override def createPubSubConnection(): GooglePubSubDAO = new SuccessfulMockGooglePubSubDao
  }

  /** A variant of PubSubMetadataServiceActor with a GooglePubSubDAO which will always return failure */
  class FailingToCreateTopicMockPubSubMetadataServiceActor(serviceConfig: Config, globalConfig: Config)
    extends PubSubMetadataServiceActor(serviceConfig, globalConfig) {

    override def createPubSubConnection(): GooglePubSubDAO = new FailingToCreateTopicMockGooglePubSubDao
  }

  /** A variant of PubSubMetadataServiceActor which will fail on message publication */
  class FailToPublishMockPubSubMetadataServiceActor(serviceConfig: Config, globalConfig: Config)
    extends PubSubMetadataServiceActor(serviceConfig, globalConfig) {

    override def createPubSubConnection(): GooglePubSubDAO = new FailToPublishMockGooglePubSubDao
  }

  trait MockGooglePubSubDao extends GooglePubSubDAO {
    override implicit val executionContext = ExecutionContext.global

    override def createTopic(topicName: String): Future[Boolean]
    override def createSubscription(topicName: String, subscriptionName: String): Future[Boolean]
    override def publishMessages(topicName: String, messages: Seq[String]): Future[Unit]

    // The following aren't used so leaving them empty
    override def deleteTopic(topicName: String): Future[Boolean] = ???
    override def getTopic(topicName: String)(implicit executionContext: ExecutionContext): Future[Option[Topic]] = ???
    override def deleteSubscription(subscriptionName: String): Future[Boolean] = ???
    override def acknowledgeMessages(subscriptionName: String, messages: Seq[PubSubMessage]): Future[Unit] = ???
    override def acknowledgeMessagesById(subscriptionName: String, ackIds: Seq[String]): Future[Unit] = ???
    override def pullMessages(subscriptionName: String, maxMessages: Int): Future[Seq[PubSubMessage]] = ???
    override def getPubSubServiceAccountCredential: Credential = ???
  }

  class SuccessfulMockGooglePubSubDao extends MockGooglePubSubDao {
    override def createTopic(topicName: String): Future[Boolean] = Future.successful(true)
    override def createSubscription(topicName: String, subscriptionName: String): Future[Boolean] = Future.successful(true)
    override def publishMessages(topicName: String, messages: Seq[String]): Future[Unit] = Future.successful(())
  }

  class FailingToCreateTopicMockGooglePubSubDao extends SuccessfulMockGooglePubSubDao {
    override def createTopic(topicName: String): Future[Boolean] = Future.failed(new RuntimeException("Unable to create topic"))
  }

  class FailToPublishMockGooglePubSubDao extends SuccessfulMockGooglePubSubDao {
    override def publishMessages(topicName: String, messages: Seq[String]): Future[Unit] = Future.failed(new RuntimeException("sorry charlie"))
  }

  // This doesn't include a project so should be a failure
  val emptyConfig = ConfigFactory.parseString(
    """
      |
    """.stripMargin
  )

  val configWithSubscription = ConfigFactory.parseString(
    """
      |project = "foo"
      |topic = "bar"
      |subscription = "baz"
    """.stripMargin
  )

  val configWithoutSubscription = ConfigFactory.parseString(
    """
      |project = "foo"
      |topic = "bar"
    """.stripMargin
  )

  val Event = MetadataEvent(MetadataKey(WorkflowId.randomId(), None, "key"), Option(MetadataValue("value")), OffsetDateTime.now)
}
