package cromwell.services


import java.util.UUID

import akka.actor.SupervisorStrategy.Stop
import akka.actor.{Actor, ActorInitializationException, ActorRef, OneForOneStrategy, Props, SupervisorStrategy}
import akka.routing.Listen
import akka.testkit.TestProbe
import com.typesafe.config.{Config, ConfigException, ConfigFactory}
import cromwell.core.TestKitSuite
import cromwell.services.BarServiceActor.{ArbitraryBarMessage, ListenToBarMessage, SetProbe}
import cromwell.services.ServiceRegistryActor.{ListenToMessage, ServiceRegistryFailure, ServiceRegistryMessage}
import cromwell.services.ServiceRegistryActorSpec._
import org.scalatest.BeforeAndAfterAll
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import scala.concurrent.duration._
import scala.language.postfixOps


abstract class EmptyActor extends Actor {
  override def receive: Receive = Actor.emptyBehavior
}

//noinspection ScalaUnusedSymbol
class FooServiceActor(config: Config, configp: Config, registryActor: ActorRef) extends EmptyActor

class NoAppropriateConstructorServiceActor extends EmptyActor

object BarServiceActor {
  trait BarServiceMessage extends ServiceRegistryMessage {
    override def serviceName = "Bar"
  }
  final case class SetProbe(probe: TestProbe) extends BarServiceMessage
  case object ArbitraryBarMessage extends BarServiceMessage
  case object ListenToBarMessage extends BarServiceMessage with ListenToMessage
}

//noinspection ScalaUnusedSymbol
class BarServiceActor(config: Config, configp: Config, registryActor: ActorRef) extends Actor {
  var probe: TestProbe = _
  override def receive: Receive = {
    case SetProbe(p) => probe = p
    case m => probe.ref forward m
  }
}

object ServiceRegistryActorSpec {

  val ServicesBlockKey = "[SERVICES_BLOCK]"

  val ConfigurationTemplate: String =
    """
      |[SERVICES_BLOCK]
    """.stripMargin

  private val ServiceNameKey = "[SERVICE_NAME]"
  private val ServiceClassKey = "[SERVICE_CLASS_NAME]"

  implicit class EnhancedServiceClass(val serviceClass: Class[_]) extends AnyVal {
    def serviceName: String = serviceClass.getSimpleName.replace("ServiceActor", "")

    def configEntry: String =
      s"""
         | $ServiceNameKey {
         |   class = "$ServiceClassKey"
         | }
      """.stripMargin.replace(ServiceNameKey, serviceClass.serviceName)
        .stripMargin.replace(ServiceClassKey, serviceClass.getCanonicalName)
  }

  val AwaitTimeout: FiniteDuration = 5 seconds
}

class ServiceRegistryActorSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with BeforeAndAfterAll {

  private def buildConfig(serviceClass: Class[_]): String = {
    val serviceEntriesKey = "[SERVICE_ENTRIES]"
    val template =
      s"""
         |  services {
         |    $serviceEntriesKey
         |  }
         |
      """.stripMargin

    val servicesBlock = template.replace(serviceEntriesKey, serviceClass.configEntry)
    ConfigurationTemplate.replace(ServicesBlockKey, servicesBlock)
  }

  private def buildProbeForInitializationException(config: Config): TestProbe = {
    val uuid = UUID.randomUUID()
    val parentProbe = TestProbe(s"parentProbe-$uuid")
    // c/p http://stackoverflow.com/questions/18619691/failing-a-scalatest-when-akka-actor-throws-exception-outside-of-the-test-thread/21892677
    system.actorOf(
      props = Props(new EmptyActor {
        context.actorOf(ServiceRegistryActor.props(config), s"ServiceRegistryActor-$uuid")
        override val supervisorStrategy: SupervisorStrategy = OneForOneStrategy() {
          case f => parentProbe.ref ! f; Stop
        }
      }),
      name = s"childActor-$uuid",
    )
    parentProbe
  }

  private def buildServiceRegistry(config: Config): ActorRef = {
    system.actorOf(ServiceRegistryActor.props(config), s"ServiceRegistryActor-${UUID.randomUUID()}")
  }

  behavior of "ServiceRegistryActorSpec"

  // The "die during construction" tests assert that an `ActorInitializationException` is sent to the Actor's supervisor
  // in a variety of circumstances which should probably trigger a Cromwell server shutdown.
  // These are unit tests and not integration tests, so there are no assertions here that the ServiceRegistryActor's
  // supervisor in the Cromwell actor hierarchy does anything sensible based upon receiving this exception.

  it should "die during construction if a config lacks a services block" in {
    val configString = buildConfig(classOf[FooServiceActor])
    val missingServices = configString.replace("  services ", "  shmervices")
    val probe = buildProbeForInitializationException(ConfigFactory.parseString(missingServices))
    probe.expectMsgPF(AwaitTimeout) {
      case e: ActorInitializationException =>
        e.getCause shouldBe a [ConfigException.Missing]
        e.getCause.getMessage shouldBe "String: 1: No configuration setting found for key 'services'"
    }
  }

  it should "die during construction if a service class can't be found" in {
    val configString = buildConfig(classOf[FooServiceActor])
    val missingService = configString.replace("FooServiceActor", "FooWhoServiceActor")
    val probe = buildProbeForInitializationException(ConfigFactory.parseString(missingService))
    probe.expectMsgPF(AwaitTimeout) {
      case e: ActorInitializationException =>
        // The class not found exception is wrapped in a Runtime Exception giving the name of the faulty service
        val cause = e.getCause
        cause shouldBe a [RuntimeException] 
        val classNotFound = cause.getCause
        classNotFound shouldBe a [ClassNotFoundException]
        classNotFound.getMessage shouldBe "cromwell.services.FooWhoServiceActor"
    }
  }

  it should "die during construction if a service class lacks a proper constructor" in {
    val configString = buildConfig(classOf[NoAppropriateConstructorServiceActor])
    val probe = buildProbeForInitializationException(ConfigFactory.parseString(configString))
    probe.expectMsgPF(AwaitTimeout) {
      case e: ActorInitializationException =>
        e.getCause shouldBe an [IllegalArgumentException]
        e.getCause.getMessage should include("no matching constructor found on class cromwell.services.NoAppropriateConstructorServiceActor")
    }
  }

  it should "die during construction if a service configuration has no 'class' attribute" in {
    val configString = buildConfig(classOf[FooServiceActor])
    val missingService = configString.replace("class = \"cromwell.services.FooServiceActor\"", "")
    val probe = buildProbeForInitializationException(ConfigFactory.parseString(missingService))
    probe.expectMsgPF(AwaitTimeout) {
      case e: ActorInitializationException =>
        e.getCause shouldBe an [IllegalArgumentException]
        e.getCause.getMessage shouldBe "Invalid configuration for service Foo: missing 'class' definition"
    }
  }

  it should "respond with a failure for a message which is not a ServiceRegistryMessage" in {
    val configString = buildConfig(classOf[FooServiceActor])
    val service = buildServiceRegistry(ConfigFactory.parseString(configString))
    val probe = TestProbe("probe")

    service.tell("This is a String, not an appropriate ServiceRegistryActor message", probe.ref)

    probe.expectMsgPF(AwaitTimeout) {
      case e: ServiceRegistryFailure =>
        e.serviceName should include("Message is not a ServiceRegistryMessage:")
    }
  }

  it should "reject a message for an unregistered service" in {
    // Configure only for the "Foo" service, send a "Bar" message.
    val configString = buildConfig(classOf[FooServiceActor])
    val service = buildServiceRegistry(ConfigFactory.parseString(configString))
    val probe = TestProbe("probe")

    service.tell(ArbitraryBarMessage, probe.ref)

    probe.expectMsgPF(AwaitTimeout) {
      case e: ServiceRegistryFailure =>
        e.serviceName shouldBe "Bar"
    }
  }

  it should "forward a message to the appropriate service" in {
    val configString = buildConfig(classOf[BarServiceActor])
    val service = buildServiceRegistry(ConfigFactory.parseString(configString))

    val barProbe = TestProbe("barProbe")
    service ! SetProbe(barProbe)
    service ! ArbitraryBarMessage

    barProbe.expectMsgPF(AwaitTimeout) {
      case ArbitraryBarMessage =>
      case x => fail("unexpected message: " + x)
    }
  }
  
  it should "unpack ListenTo messages" in {
    val snd = TestProbe("snd").ref
    val configString = buildConfig(classOf[BarServiceActor])
    val service = buildServiceRegistry(ConfigFactory.parseString(configString))

    val barProbe = TestProbe("barProbe")
    service ! SetProbe(barProbe)
    service.tell(ListenToBarMessage, snd)

    barProbe.expectMsgPF(AwaitTimeout) {
      case Listen(l) => l shouldBe snd
      case x => fail("unexpected message: " + x)
    }
  }
}
