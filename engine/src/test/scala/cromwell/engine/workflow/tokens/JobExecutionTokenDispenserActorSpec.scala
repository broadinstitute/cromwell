package cromwell.engine.workflow.tokens

import java.util.UUID

import akka.actor.{ActorSystem, PoisonPill}
import akka.testkit.{ImplicitSender, TestActorRef, TestKit, TestProbe}
import cromwell.core.JobExecutionToken
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor.{JobExecutionTokenDenied, JobExecutionTokenDispensed, JobExecutionTokenRequest, JobExecutionTokenReturn}
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActorSpec._
import cromwell.engine.workflow.tokens.TestTokenGrabbingActor.StoppingSupervisor
import cromwell.util.AkkaTestUtil
import org.scalatest._
import org.scalatest.concurrent.Eventually

import scala.concurrent.duration._

class JobExecutionTokenDispenserActorSpec extends TestKit(ActorSystem("JETDASpec")) with ImplicitSender with FlatSpecLike with Matchers with BeforeAndAfter with BeforeAndAfterAll with Eventually {

  val MaxWaitTime = 100.milliseconds
  implicit val pc: PatienceConfig = PatienceConfig(MaxWaitTime)

  behavior of "JobExecutionTokenDispenserActor"

  it should "dispense an infinite token correctly" in {
    actorRefUnderTest ! JobExecutionTokenRequest(TestInfiniteTokenType)
    expectMsgPF(max = MaxWaitTime, hint = "token dispensed message") {
      case JobExecutionTokenDispensed(token) =>
        token.jobExecutionTokenType should be(TestInfiniteTokenType)
    }
  }

  it should "accept return of an infinite token correctly" in {
    actorRefUnderTest ! JobExecutionTokenRequest(TestInfiniteTokenType)
    expectMsgPF(max = MaxWaitTime, hint = "token dispensed message") {
      case JobExecutionTokenDispensed(token) =>
        actorRefUnderTest ! JobExecutionTokenReturn(token)
    }
  }

  it should "dispense indefinitely for an infinite token type" in {
    var currentSet: Set[JobExecutionToken] = Set.empty
    100 indexedTimes { i =>
      val sender = TestProbe()
      actorRefUnderTest.tell(msg = JobExecutionTokenRequest(TestInfiniteTokenType), sender = sender.ref)
      sender.expectMsgPF(max = MaxWaitTime, hint = "token dispensed message") {
        case JobExecutionTokenDispensed(token) =>
          token.jobExecutionTokenType should be(TestInfiniteTokenType)
          currentSet.contains(token) should be(false)
          currentSet += token
      }
    }
  }

  it should "dispense a limited token correctly" in {

    actorRefUnderTest ! JobExecutionTokenRequest(LimitedTo5Tokens)
    expectMsgPF(max = MaxWaitTime, hint = "token dispensed message") {
      case JobExecutionTokenDispensed(token) => token.jobExecutionTokenType should be(LimitedTo5Tokens)
    }
  }

  it should "accept return of a limited token type correctly" in {
    actorRefUnderTest ! JobExecutionTokenRequest(LimitedTo5Tokens)
    expectMsgPF(max = MaxWaitTime, hint = "token dispensed message") {
      case JobExecutionTokenDispensed(token) => actorRefUnderTest ! JobExecutionTokenReturn(token)
    }
  }

  it should "limit the dispensing of a limited token type" in {

    var currentTokens: Map[TestProbe, JobExecutionToken] = Map.empty
    val dummyActors = (0 until 100 map { i => i -> TestProbe("dummy_" + i) }).toMap

    // Dispense the first 5:
    5 indexedTimes { i =>
      val sndr = dummyActors(i)
      actorRefUnderTest.tell(msg = JobExecutionTokenRequest(LimitedTo5Tokens), sender = sndr.ref)
      sndr.expectMsgPF(max = MaxWaitTime, hint = "token dispensed message") {
        case JobExecutionTokenDispensed(token) =>
          token.jobExecutionTokenType should be(LimitedTo5Tokens)
          currentTokens.values.toList.contains(token) should be(false) // Check we didn't already get this token
          currentTokens += sndr -> token
      }
    }

    // Queue the next 95:
    95 indexedTimes { i =>
      val sndr = dummyActors(5 + i)
      actorRefUnderTest.tell(msg = JobExecutionTokenRequest(LimitedTo5Tokens), sender = sndr.ref)
      sndr.expectMsgPF(max = MaxWaitTime, hint = "token denied message") {
        case JobExecutionTokenDenied(positionInQueue) =>
          positionInQueue should be(i)
      }
    }

    // It should allow queued actors to check their position in the queue:
    95 indexedTimes { i =>
      val sndr = dummyActors(5 + i)
      actorRefUnderTest.tell(msg = JobExecutionTokenRequest(LimitedTo5Tokens), sender = sndr.ref)
      sndr.expectMsgPF(max = MaxWaitTime, hint = "token denied message") {
        case JobExecutionTokenDenied(positionInQueue) =>
          positionInQueue should be(i)
      }
    }

    // It should release tokens as soon as they're available (while there's still a queue...):
    95 indexedTimes { i =>
      val returner = dummyActors(i)
      val nextInLine = dummyActors(i + 5)
      val tokenBeingReturned = currentTokens(returner)
      actorRefUnderTest.tell(msg = JobExecutionTokenReturn(tokenBeingReturned), sender = returner.ref)
      currentTokens -= returner
      nextInLine.expectMsgPF(max = MaxWaitTime, hint = s"token dispensed message to the next in line actor (#${i + 5})") {
        case JobExecutionTokenDispensed(token) =>
          token should be(tokenBeingReturned) // It just gets immediately passed out again!
          currentTokens += nextInLine -> token
      }
    }

    // Double-check the queue state: when we request a token now, we should still be denied:
    actorRefUnderTest ! JobExecutionTokenRequest(LimitedTo5Tokens)
    expectMsgClass(classOf[JobExecutionTokenDenied])

    //And finally, silently release the remaining tokens:
    5 indexedTimes { i =>
      val returner = dummyActors(i + 95)
      val tokenBeingReturned = currentTokens(returner)
      actorRefUnderTest.tell(msg = JobExecutionTokenReturn(tokenBeingReturned), sender = returner.ref)
      currentTokens -= returner
    }

    // And we should have gotten our own token by now:
    expectMsgClass(classOf[JobExecutionTokenDispensed])

    // Check we didn't get anything else in the meanwhile:
    msgAvailable should be(false)
    dummyActors.values foreach { testProbe => testProbe.msgAvailable should be(false) }
  }

  it should "resend the same token to an actor which already has one" in {
    actorRefUnderTest ! JobExecutionTokenRequest(LimitedTo5Tokens)
    val firstResponse = expectMsgClass(classOf[JobExecutionTokenDispensed])

    5 indexedTimes { i =>
      actorRefUnderTest ! JobExecutionTokenRequest(LimitedTo5Tokens)
      expectMsg(MaxWaitTime, s"same token again (attempt ${i + 1})", firstResponse) // Always the same
    }
  }


  // Incidentally, also covers: it should "not be fooled if the wrong actor returns a token"
  it should "not be fooled by a doubly-returned token" in {
    var currentTokens: Map[TestProbe, JobExecutionToken] = Map.empty
    val dummyActors = (0 until 7 map { i => i -> TestProbe("dummy_" + i) }).toMap

    // Set up by taking all 5 tokens out, and then adding 2 to the queue:
    5 indexedTimes { i =>
      val sndr = dummyActors(i)
      actorRefUnderTest.tell(msg = JobExecutionTokenRequest(LimitedTo5Tokens), sender = sndr.ref)
      currentTokens += dummyActors(i) -> sndr.expectMsgClass(classOf[JobExecutionTokenDispensed]).jobExecutionToken
    }
    2 indexedTimes { i =>
      val sndr = dummyActors(5 + i)
      actorRefUnderTest.tell(msg = JobExecutionTokenRequest(LimitedTo5Tokens), sender = sndr.ref)
      sndr.expectMsgClass(classOf[JobExecutionTokenDenied])
    }

    // The first time we return a token, the next in line should be given it:
    val returningActor = dummyActors(0)
    val nextInLine1 = dummyActors(5)
    val nextInLine2 = dummyActors(6)
    val tokenBeingReturned = currentTokens(returningActor)
    currentTokens -= returningActor
    actorRefUnderTest.tell(msg = JobExecutionTokenReturn(tokenBeingReturned), sender = returningActor.ref)
    val tokenPassedOn = nextInLine1.expectMsgClass(classOf[JobExecutionTokenDispensed]).jobExecutionToken
    tokenPassedOn should be(tokenBeingReturned)
    currentTokens += nextInLine1 -> tokenPassedOn

    // But the next time, nothing should happen because the wrong actor is returning the token:
    actorRefUnderTest.tell(msg = JobExecutionTokenReturn(tokenBeingReturned), sender = returningActor.ref)
    nextInLine2.expectNoMsg(MaxWaitTime)
  }

  it should "not be fooled if an actor returns a token which doesn't exist" in {
    var currentTokens: Map[TestProbe, JobExecutionToken] = Map.empty
    val dummyActors = (0 until 6 map { i => i -> TestProbe("dummy_" + i) }).toMap

    // Set up by taking all 5 tokens out, and then adding 2 to the queue:
    5 indexedTimes { i =>
      val sndr = dummyActors(i)
      actorRefUnderTest.tell(msg = JobExecutionTokenRequest(LimitedTo5Tokens), sender = sndr.ref)
      currentTokens += dummyActors(i) -> sndr.expectMsgClass(classOf[JobExecutionTokenDispensed]).jobExecutionToken
    }
    1 indexedTimes { i =>
      val sndr = dummyActors(5 + i)
      actorRefUnderTest.tell(msg = JobExecutionTokenRequest(LimitedTo5Tokens), sender = sndr.ref)
      sndr.expectMsgClass(classOf[JobExecutionTokenDenied])
    }

    actorRefUnderTest.tell(msg = JobExecutionTokenReturn(JobExecutionToken(LimitedTo5Tokens, UUID.randomUUID())), sender = dummyActors(0).ref)
    dummyActors(5).expectNoMsg(MaxWaitTime)
  }

  AkkaTestUtil.actorDeathMethods(system) foreach { case (name, stopMethod) =>
    it should s"recover tokens lost to actors which are $name before they hand back their token" in {
      var currentTokens: Map[TestActorRef[TestTokenGrabbingActor], JobExecutionToken] = Map.empty
      var tokenGrabbingActors: Map[Int, TestActorRef[TestTokenGrabbingActor]] = Map.empty
      val grabberSupervisor = TestActorRef(new StoppingSupervisor())

      // Set up by taking all 5 tokens out, and then adding 2 to the queue:
      5 indexedTimes { i =>
        val newGrabbingActor = TestActorRef[TestTokenGrabbingActor](TestTokenGrabbingActor.props(actorRefUnderTest, LimitedTo5Tokens), grabberSupervisor, s"grabber_${name}_" + i)
        tokenGrabbingActors += i -> newGrabbingActor
        eventually {
          newGrabbingActor.underlyingActor.token.isDefined should be(true)
        }
        currentTokens += newGrabbingActor -> newGrabbingActor.underlyingActor.token.get
      }

      val unassignedActorIndex = 5
      val newGrabbingActor = TestActorRef(new TestTokenGrabbingActor(actorRefUnderTest, LimitedTo5Tokens), s"grabber_${name}_" + unassignedActorIndex)
      tokenGrabbingActors += unassignedActorIndex -> newGrabbingActor
      eventually {
        newGrabbingActor.underlyingActor.rejections should be(1)
      }

      val actorToStop = tokenGrabbingActors(0)
      val actorToStopsToken = currentTokens(actorToStop)
      val nextInLine = tokenGrabbingActors(unassignedActorIndex)

      val deathwatch = TestProbe()
      deathwatch watch actorToStop
      stopMethod(actorToStop)
      deathwatch.expectTerminated(actorToStop)
      eventually { nextInLine.underlyingActor.token should be(Some(actorToStopsToken)) }
    }
  }

  it should "skip over dead actors when assigning tokens to the actor queue" in {
    var currentTokens: Map[TestActorRef[TestTokenGrabbingActor], JobExecutionToken] = Map.empty
    var tokenGrabbingActors: Map[Int, TestActorRef[TestTokenGrabbingActor]] = Map.empty
    val grabberSupervisor = TestActorRef(new StoppingSupervisor())

    // Set up by taking all 5 tokens out, and then adding 2 to the queue:
    5 indexedTimes { i =>
      val newGrabbingActor = TestActorRef[TestTokenGrabbingActor](TestTokenGrabbingActor.props(actorRefUnderTest, LimitedTo5Tokens), grabberSupervisor, s"skip_test_" + i)
      tokenGrabbingActors += i -> newGrabbingActor
      eventually {
        newGrabbingActor.underlyingActor.token.isDefined should be(true)
      }
      currentTokens += newGrabbingActor -> newGrabbingActor.underlyingActor.token.get
    }
    2 indexedTimes { i =>
      val index = i + 5
      val newGrabbingActor = TestActorRef[TestTokenGrabbingActor](TestTokenGrabbingActor.props(actorRefUnderTest, LimitedTo5Tokens), grabberSupervisor, s"skip_test_" + index)
      tokenGrabbingActors += index -> newGrabbingActor
      eventually {
        newGrabbingActor.underlyingActor.rejections should be(1)
      }
    }

    val returningActor = tokenGrabbingActors(0)
    val returnedToken = currentTokens(returningActor)
    val nextInLine1 = tokenGrabbingActors(5)
    val nextInLine2 = tokenGrabbingActors(6)

    // First, kill off the actor which would otherwise be first in line:
    val deathwatch = TestProbe()
    deathwatch watch nextInLine1
    nextInLine1 ! PoisonPill
    deathwatch.expectTerminated(nextInLine1)

    // Now, stop one of the workers unexpectedly and check that the released token goes to the right place:
    actorRefUnderTest.tell(msg = JobExecutionTokenReturn(returnedToken), sender = returningActor)
    eventually { nextInLine2.underlyingActor.token should be(Some(returnedToken)) } // Some is OK. This is the **expected** value!
  }

  var actorRefUnderTest: TestActorRef[JobExecutionTokenDispenserActor] = _

  before {
    actorRefUnderTest = TestActorRef(new JobExecutionTokenDispenserActor())

  }
  after {
    actorRefUnderTest = null
  }

  override def afterAll = {
    TestKit.shutdownActorSystem(system)
  }
}

object JobExecutionTokenDispenserActorSpec {

  implicit class intWithTimes(n: Int) {
    def times(f: => Unit) = 1 to n foreach { _ => f }
    def indexedTimes(f: Int => Any) = 0 until n foreach { i => f(i) }
  }

  val TestInfiniteTokenType = JobExecutionTokenType("infinite", maxPoolSize = None)
  def limitedTokenType(limit: Int) = JobExecutionTokenType(s"$limit-limit", maxPoolSize = Option(limit))
  val LimitedTo5Tokens = limitedTokenType(5)
}
