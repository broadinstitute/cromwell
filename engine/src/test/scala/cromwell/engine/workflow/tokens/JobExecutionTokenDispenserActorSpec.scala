package cromwell.engine.workflow.tokens

import akka.actor.{ActorRef, ActorSystem, PoisonPill}
import akka.testkit.{ImplicitSender, TestActorRef, TestKit, TestProbe}
import cromwell.core.HogGroup
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.DynamicRateLimiter.{Rate, TokensAvailable}
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor._
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActorSpec._
import cromwell.engine.workflow.tokens.TokenQueue.TokenQueuePlaceholder
import cromwell.util.AkkaTestUtil
import cromwell.util.AkkaTestUtil.StoppingSupervisor
import org.scalatest._
import org.scalatest.concurrent.Eventually

import scala.concurrent.duration._
import scala.util.Random

class JobExecutionTokenDispenserActorSpec extends TestKit(ActorSystem("JETDASpec")) with ImplicitSender with FlatSpecLike with Matchers with BeforeAndAfter with BeforeAndAfterAll with Eventually {

  val MaxWaitTime = 10.seconds
  implicit val pc: PatienceConfig = PatienceConfig(MaxWaitTime)

  behavior of "JobExecutionTokenDispenserActor"

  val hogGroupA = HogGroup("hogGroupA")
  val hogGroupB = HogGroup("hogGroupB")

  it should "dispense an infinite token correctly" in {
    actorRefUnderTest ! JobExecutionTokenRequest(hogGroupA, TestInfiniteTokenType)
    expectMsg(max = MaxWaitTime, JobExecutionTokenDispensed)
    actorRefUnderTest.underlyingActor.tokenAssignments.contains(self) shouldBe true
    actorRefUnderTest.underlyingActor.tokenAssignments(self).get().jobExecutionTokenType shouldBe TestInfiniteTokenType
  }

  it should "accept return of an infinite token correctly" in {
    actorRefUnderTest ! JobExecutionTokenRequest(hogGroupA, TestInfiniteTokenType)
    expectMsg(max = MaxWaitTime, JobExecutionTokenDispensed)
    actorRefUnderTest.underlyingActor.tokenAssignments.contains(self) shouldBe true
    actorRefUnderTest.underlyingActor.tokenAssignments(self).get().jobExecutionTokenType shouldBe TestInfiniteTokenType
    actorRefUnderTest ! JobExecutionTokenReturn
    actorRefUnderTest.underlyingActor.tokenAssignments.contains(self) shouldBe false
  }

  it should "dispense indefinitely for an infinite token type" in {
    val senders = (1 to 20).map(_ => TestProbe())
    senders.foreach(sender => actorRefUnderTest.tell(msg = JobExecutionTokenRequest(hogGroupA, TestInfiniteTokenType), sender = sender.ref))
    senders.foreach(_.expectMsg(max = MaxWaitTime, JobExecutionTokenDispensed))
    actorRefUnderTest.underlyingActor.tokenAssignments.size shouldBe 20
  }

  it should "dispense the correct amount at the specified rate, not more and not faster" in {

    // Override with a slower distribution rate for this one test:
    actorRefUnderTest = TestActorRef(new JobExecutionTokenDispenserActor(TestProbe().ref, Rate(10, 4.seconds), None))

    val senders = (1 to 20).map(_ => TestProbe())
    senders.foreach(sender => actorRefUnderTest.tell(msg = JobExecutionTokenRequest(hogGroupA, TestInfiniteTokenType), sender = sender.ref))
    // The first 10 should get their token
    senders.take(10).foreach(_.expectMsg(max = MaxWaitTime, JobExecutionTokenDispensed))
    // Couldn't figure out a cleaner way to "verify that none of this probes gets a message in the next X seconds"
    Thread.sleep(1.second.toMillis)
    // Then the last ten should eventually get it too, but later
    senders.drop(10).foreach(_.msgAvailable shouldBe false)
    senders.slice(10, 20).foreach(_.expectMsg(max = MaxWaitTime, JobExecutionTokenDispensed))
  }

  it should "dispense a limited token correctly" in {
    actorRefUnderTest ! JobExecutionTokenRequest(hogGroupA, LimitedTo5Tokens)
    expectMsg(max = MaxWaitTime, JobExecutionTokenDispensed)
    actorRefUnderTest.underlyingActor.tokenAssignments.contains(self) shouldBe true
    actorRefUnderTest.underlyingActor.tokenAssignments(self).get().jobExecutionTokenType shouldBe LimitedTo5Tokens
  }

  it should "accept return of a limited token type correctly" in {
    actorRefUnderTest ! JobExecutionTokenRequest(hogGroupA, LimitedTo5Tokens)
    expectMsg(max = MaxWaitTime, JobExecutionTokenDispensed)
    actorRefUnderTest.underlyingActor.tokenAssignments.contains(self) shouldBe true
    actorRefUnderTest.underlyingActor.tokenAssignments(self).get().jobExecutionTokenType shouldBe LimitedTo5Tokens
    actorRefUnderTest ! JobExecutionTokenReturn
    actorRefUnderTest.underlyingActor.tokenAssignments.contains(self) shouldBe false
  }

  def asHogGroupAPlaceholder(probe: TestProbe): TokenQueuePlaceholder = asHogGroupAPlaceholder(probe.ref)
  def asHogGroupAPlaceholder(ref: ActorRef): TokenQueuePlaceholder = TokenQueuePlaceholder(ref, "hogGroupA")

  it should "limit the dispensing of a limited token type" in {
    val senders = (1 to 15).map(_ => TestProbe())
    // Ask for 20 tokens
    senders.foreach(sender => actorRefUnderTest.tell(msg = JobExecutionTokenRequest(hogGroupA, LimitedTo5Tokens), sender = sender.ref))

    // Force token distribution
    actorRefUnderTest ! TokensAvailable(100)
    senders.take(5).foreach(_.expectMsg(JobExecutionTokenDispensed))

    actorRefUnderTest.underlyingActor.tokenAssignments.size shouldBe 5
    actorRefUnderTest.underlyingActor.tokenAssignments.keySet should contain theSameElementsAs senders.map(_.ref).take(5).toSet

    // The last 10 should be queued
    // At this point [0, 1, 2, 3, 4] are the ones with tokens, and [5, 14] are still queued
    actorRefUnderTest.underlyingActor.tokenQueues(LimitedTo5Tokens).size shouldBe 10
    actorRefUnderTest.underlyingActor.tokenQueues(LimitedTo5Tokens).queues.flatMap(_._2).toList should contain theSameElementsInOrderAs senders.drop(5).map(asHogGroupAPlaceholder)

    // Force token distribution
    actorRefUnderTest ! TokensAvailable(100)
    // The other still should have received nothing
    senders.drop(5).foreach(_.msgAvailable shouldBe false)

    // Release a few tokens
    // The next 3 (i.e 5, 6, 7) should get their token at the next distribution cycle
    senders.take(3).foreach(_.send(actorRefUnderTest, JobExecutionTokenReturn))

    // Force token distribution
    actorRefUnderTest ! TokensAvailable(100)
    senders.slice(5, 8).foreach(_.expectMsg(JobExecutionTokenDispensed))

    // At this point [3, 4, 5, 6, 7] are the ones with tokens, and [8, 19] are still queued
    actorRefUnderTest.underlyingActor.tokenQueues(LimitedTo5Tokens).queues.flatMap(_._2).toList should contain theSameElementsInOrderAs senders.slice(8, 20).map(asHogGroupAPlaceholder)

    // Double-check the queue state: when we request a token now, we should still be denied:
    actorRefUnderTest ! JobExecutionTokenRequest(hogGroupA, LimitedTo5Tokens)
    // Force token distribution
    actorRefUnderTest ! TokensAvailable(100)
    expectNoMessage()
    // We should be enqueued and the last in the queue though
    actorRefUnderTest.underlyingActor.tokenQueues(LimitedTo5Tokens).queues.flatMap(_._2).last shouldBe TokenQueuePlaceholder(self, "hogGroupA")

    // Release all currently owned tokens
    senders.slice(3, 8).foreach(_.send(actorRefUnderTest, JobExecutionTokenReturn))
    // Force token distribution
    actorRefUnderTest ! TokensAvailable(100)
    actorRefUnderTest.underlyingActor.tokenAssignments.keySet should contain theSameElementsAs senders.map(_.ref).slice(8, 13)
    // Keep accepting and returning tokens immediately
    senders.slice(8, 13).foreach(_.expectMsg(JobExecutionTokenDispensed))
    senders.slice(8, 13).foreach(_.reply(JobExecutionTokenReturn))
    actorRefUnderTest ! TokensAvailable(100)
    // Last 2 tokens
    senders.slice(13, 15).foreach(_.expectMsg(JobExecutionTokenDispensed))

    // We were last on the list but we should have our token now
    expectMsg(JobExecutionTokenDispensed)

    // Queue should be empty now
    actorRefUnderTest.underlyingActor.tokenQueues(LimitedTo5Tokens).size shouldBe 0

    // There should be 3 assigned tokens: index 18, 19, and this test actor
    actorRefUnderTest.underlyingActor.tokenAssignments.keySet should contain theSameElementsAs senders.map(_.ref).slice(13, 15) :+ self
  }

  it should "resend the same token to an actor which already has one" in {
    5 indexedTimes { _ =>
      actorRefUnderTest ! JobExecutionTokenRequest(hogGroupA, LimitedTo5Tokens)
    }
    // Force token distribution
    actorRefUnderTest ! TokensAvailable(100)
    // Nothing should be queued
    actorRefUnderTest.underlyingActor.tokenQueues(LimitedTo5Tokens).size shouldBe 0
    // We should be the only actor with an assigned token
    actorRefUnderTest.underlyingActor.tokenAssignments.keySet shouldBe Set(self)
    // Only 1 token should be leased
    actorRefUnderTest.underlyingActor.tokenQueues(LimitedTo5Tokens).pool.leased() shouldBe 1
  }


  //Incidentally, also covers: it should "not be fooled if the wrong actor returns a token"
  it should "not be fooled by a doubly-returned token" in {
    val senders = (1 to 7).map(_ => TestProbe())
    // Ask for 7 tokens
    senders.foreach(sender => actorRefUnderTest.tell(msg = JobExecutionTokenRequest(hogGroupA, LimitedTo5Tokens), sender = sender.ref))

    // Force token distribution
    actorRefUnderTest ! TokensAvailable(5)
    // Get the first 5
    senders.take(5).foreach(_.expectMsg(JobExecutionTokenDispensed))

    // Sender 0 returns his token
    actorRefUnderTest.tell(JobExecutionTokenReturn, senders.head.ref)

    // Force token distribution
    actorRefUnderTest ! TokensAvailable(1)
    // Sender 5 gets his
    senders(5).expectMsg(JobExecutionTokenDispensed)

    // Now sender 0 again returns a token, although it doesn't have one !
    actorRefUnderTest.tell(JobExecutionTokenReturn, senders.head.ref)

    // Nothing should happen, specifically Sender 6 should not get his token
    senders(6).expectNoMessage(MaxWaitTime)
  }

  AkkaTestUtil.actorDeathMethods(system) foreach { case (name, stopMethod) =>
    it should s"recover tokens lost to actors which are $name before they hand back their token" in {
      val grabberSupervisor = TestActorRef(new StoppingSupervisor())
      // The first 5 get a token and the 6th one is queued
      val tokenGrabbingActors = (1 to 6).map { i =>
        TestActorRef[TestTokenGrabbingActor](TestTokenGrabbingActor.props(actorRefUnderTest, LimitedTo5Tokens), grabberSupervisor, s"grabber_${name}_" + i)
      }

      // Force token distribution
      actorRefUnderTest ! TokensAvailable(5)
      // Get the first 5
      tokenGrabbingActors.take(5).foreach(_.underlyingActor.hasToken shouldBe true)
      tokenGrabbingActors(5).underlyingActor.hasToken shouldBe false

      // Stop the first one
      val actorToStop = tokenGrabbingActors.head

      // Expect the last one to get his token
      val nextInLine = tokenGrabbingActors.last

      val deathwatch = TestProbe()
      deathwatch watch actorToStop
      stopMethod(actorToStop)
      deathwatch.expectTerminated(actorToStop)
      eventually { nextInLine.underlyingActor.hasToken shouldBe true }
    }
  }

  it should "skip over dead actors when assigning tokens to the actor queue" in {
    val grabberSupervisor = TestActorRef(new StoppingSupervisor())
    // The first 5 get a token and the 6th and 7h one are queued
    val tokenGrabbingActors = (1 to 7).map { i =>
      TestActorRef[TestTokenGrabbingActor](TestTokenGrabbingActor.props(actorRefUnderTest, LimitedTo5Tokens), grabberSupervisor, s"grabber_" + i)
    }

    // Force token distribution
    actorRefUnderTest ! TokensAvailable(5)
    // Get the first 5
    tokenGrabbingActors.take(5).foreach(_.underlyingActor.hasToken shouldBe true)
    val nextInLine1 = tokenGrabbingActors(5)
    val nextInLine2 = tokenGrabbingActors(6)

    // Check that the next in lines have no tokens and are indeed in the queue
    nextInLine1.underlyingActor.hasToken shouldBe false
    nextInLine2.underlyingActor.hasToken shouldBe false
    actorRefUnderTest.underlyingActor.tokenQueues(LimitedTo5Tokens).queues.flatMap(_._2).toList should contain theSameElementsInOrderAs List(nextInLine1, nextInLine2).map(asHogGroupAPlaceholder)

    // First, kill off the actor which would otherwise be first in line:
    val deathwatch = TestProbe()
    deathwatch watch nextInLine1
    nextInLine1 ! PoisonPill
    deathwatch.expectTerminated(nextInLine1)

    // Now, have an actor return its token and check that the released token goes to nextInLine2:
    actorRefUnderTest.tell(msg = JobExecutionTokenReturn, sender = tokenGrabbingActors.head)
    // Force token distribution
    actorRefUnderTest ! TokensAvailable(1)
    eventually { nextInLine2.underlyingActor.hasToken shouldBe true }
  }

  it should "skip over dead actors repeatedly when assigning tokens to the actor queue" in {
    val grabberSupervisor = TestActorRef(new StoppingSupervisor())
    // The first 5 get a token and the 6th and 7h one are queued
    val tokenGrabbingActors = (0 until 1000).toVector.map { i =>
      TestActorRef[TestTokenGrabbingActor](TestTokenGrabbingActor.props(actorRefUnderTest, LimitedTo5Tokens), grabberSupervisor, s"grabber_" + i)
    }

    val actorIterator = tokenGrabbingActors.toIterator

    while (actorIterator.hasNext) {

      // We won't actually dispense 100, this is simulating the "steady drip" message
      // so that we don't have to wait 4 seconds per drip for the test case...
      actorRefUnderTest ! TokensAvailable(100)
      val withTokens = actorIterator.take(5).toList
      val nextInLine = actorIterator.take(5).toList

      eventually {
        withTokens.foreach(_.underlyingActor.hasToken should be(true))
      }

      // Check that the next in lines have no tokens and are indeed in the queue
      nextInLine.foreach(next => next.underlyingActor.hasToken shouldBe false)

      // Kill off running jobs and the queuers (except the 4th in line):
      (0 until 5).filterNot(_ == 3) foreach { index =>
        // Kill off jobs in the queue:
        nextInLine(index) ! PoisonPill
        // Stop or kill off jobs which are 'running'
        val randomInt = Random.nextInt(10)
        if (randomInt <= 5) withTokens(index) ! PoisonPill
        else actorRefUnderTest.tell(msg = JobExecutionTokenReturn, sender = withTokens(index))
      }
      // Also complete the 4th running job (but not the 4th in the queue):
      actorRefUnderTest.tell(msg = JobExecutionTokenReturn, sender = withTokens(3))

      actorRefUnderTest ! TokensAvailable(100)
      eventually { nextInLine(3).underlyingActor.hasToken shouldBe true }

      // And kill off the rest of the actors:
      (withTokens :+ nextInLine(3)) foreach { actor => actor ! PoisonPill }
    }

    actorRefUnderTest.underlyingActor.tokenQueues.map(x => x._2.size).sum should be(0)
  }

  it should "be resilient if the last request for a hog group is removed from the queue before a token is dispensed" in {
    val tokenType = JobExecutionTokenType(s"mini", maxPoolSize = Option(6), hogFactor = 2)

    val groupARequesters = (0 until 4) map { i => TestProbe(name = s"group_a_$i") }
    val groupBRequesters = (0 until 4) map { i => TestProbe(name = s"group_b_$i") }

    // Group A and B fill up the token queue:
    groupARequesters foreach { probe => probe.send(actorRefUnderTest, JobExecutionTokenRequest(hogGroupA, tokenType)) }
    groupBRequesters foreach { probe => probe.send(actorRefUnderTest, JobExecutionTokenRequest(hogGroupB, tokenType)) }

    // Pretty soon, the tokens should all be gone, 3 to each group, and each group should have one queued item:
    eventually {
      (groupARequesters.take(3) ++ groupBRequesters.take(3)).foreach { probe =>
        actorRefUnderTest.underlyingActor.tokenAssignments.keys should contain(probe.ref)
      }
      // And both groups should have one item in the queue:
      actorRefUnderTest.underlyingActor.tokenQueues(tokenType).queues.values.foreach { queue => queue.size should be(1) }
    }

    // Group B gets bored and aborts all jobs (in reverse order to make sure ):
    groupBRequesters.reverse.foreach { probe =>
      probe.ref ! PoisonPill
      // TODO: validate that the probe has left the queue
    }

    // Group A's jobs are able to complete successfully:
    groupARequesters foreach { probe =>
      probe.expectMsg(JobExecutionTokenDispensed)
      probe.ref ! PoisonPill
    }

  }

  var actorRefUnderTest: TestActorRef[JobExecutionTokenDispenserActor] = _

  before {
    actorRefUnderTest = TestActorRef(new JobExecutionTokenDispenserActor(TestProbe().ref, Rate(10, 100.millis), None))
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

  val TestInfiniteTokenType = JobExecutionTokenType("infinite", maxPoolSize = None, hogFactor = 1)
  def limitedTokenType(limit: Int) = JobExecutionTokenType(s"$limit-limit", maxPoolSize = Option(limit), hogFactor = 1)
  val LimitedTo5Tokens = limitedTokenType(5)
}
