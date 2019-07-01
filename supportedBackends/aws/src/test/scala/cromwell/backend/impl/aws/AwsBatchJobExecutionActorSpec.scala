/*
 * Copyright 2018 Amazon.com, Inc. or its affiliates.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 *  BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 *  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 *  THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 *  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 *  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 *  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

package cromwell.backend.impl.aws

import akka.actor.{Actor, ActorRef, Props}
import akka.testkit._
import cromwell.backend.BackendJobExecutionActor.{ExecuteJobCommand, JobFailedNonRetryableResponse}
import cromwell.backend.impl.aws.ControllableFailingJabjea.JabjeaExplode
import cromwell.backend.standard.{DefaultStandardSyncExecutionActorParams, StandardSyncExecutionActor, StandardSyncExecutionActorParams}
import cromwell.backend.{BackendJobDescriptor, MinimumRuntimeSettings}
import cromwell.core.TestKitSuite
import org.scalatest.{FlatSpecLike, Matchers}
import org.specs2.mock.Mockito

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Promise}
import scala.util.control.NoStackTrace
import scala.util.{Failure, Success}

class AwsBatchJobExecutionActorSpec extends TestKitSuite("AwsBatchJobExecutionActorSpec") with FlatSpecLike with Matchers with Mockito {

  behavior of "AwsBatchJobExecutionActor"

  private val AwaitAlmostNothing = 100.milliseconds.dilated
  private val TimeoutDuration = 10.seconds.dilated
  implicit val ec: ExecutionContext = system.dispatcher

  it should "catch failures in execution actor initialization and fail the job accordingly" in {
    val jobDescriptor = BackendJobDescriptor(null, null, null, Map.empty, null, null, null)
    val workflowInfo = mock[AwsBatchConfiguration]
    val initializationData = mock[AwsBatchBackendInitializationData]
    val serviceRegistryActor = system.actorOf(Props.empty)
    val ioActor = system.actorOf(Props.empty)
    val backendSingletonActor = Option(system.actorOf(Props.empty))

    initializationData.configuration returns workflowInfo

    val parent = TestProbe()
    val deathwatch = TestProbe()
    val params = DefaultStandardSyncExecutionActorParams(AwsBatchAsyncBackendJobExecutionActor.AwsBatchOperationIdKey, serviceRegistryActor, ioActor,
      jobDescriptor, null, Option(initializationData), backendSingletonActor,
      classOf[AwsBatchAsyncBackendJobExecutionActor], MinimumRuntimeSettings())
    val testJJEA = TestActorRef[TestAwsBatchJobExecutionActor](
      props = Props(new TestAwsBatchJobExecutionActor(params, Props(new ConstructorFailingJABJEA))),
      supervisor = parent.ref)
    deathwatch watch testJJEA

    // Nothing happens:
    parent.expectNoMessage(max = AwaitAlmostNothing)
    deathwatch.expectNoMessage(max = AwaitAlmostNothing)

    testJJEA.tell(msg = ExecuteJobCommand, sender = parent.ref)

    parent.expectMsgPF(max = TimeoutDuration) {
      case JobFailedNonRetryableResponse(_, throwable, _) =>
        throwable.getMessage should be("AwsBatchAsyncBackendJobExecutionActor failed and didn't catch its exception. This condition has been handled and the job will be marked as failed.")
    }
  }

  it should "catch failures at a random point during execution actor processing and fail the job accordingly" in {
    val jobDescriptor = BackendJobDescriptor(null, null, null, Map.empty, null, null, null)
    val workflowInfo = mock[AwsBatchConfiguration]
    val initializationData = mock[AwsBatchBackendInitializationData]
    val serviceRegistryActor = system.actorOf(Props.empty)
    val ioActor = system.actorOf(Props.empty)
    val backendSingletonActor = Option(system.actorOf(Props.empty))

    initializationData.configuration returns workflowInfo

    val parent = TestProbe()
    val deathwatch = TestProbe()
    val constructionPromise = Promise[ActorRef]()
    val params = DefaultStandardSyncExecutionActorParams(AwsBatchAsyncBackendJobExecutionActor.AwsBatchOperationIdKey, serviceRegistryActor, ioActor,
      jobDescriptor, null, Option(initializationData), backendSingletonActor,
      classOf[AwsBatchAsyncBackendJobExecutionActor],
      MinimumRuntimeSettings())
    val testJJEA = TestActorRef[TestAwsBatchJobExecutionActor](
      props = Props(new TestAwsBatchJobExecutionActor(params, Props(new ControllableFailingJabjea(constructionPromise)))),
      supervisor = parent.ref)
    deathwatch watch testJJEA

    // Nothing happens:
    parent.expectNoMessage(max = AwaitAlmostNothing)
    deathwatch.expectNoMessage(max = AwaitAlmostNothing)

    testJJEA.tell(msg = ExecuteJobCommand, sender = parent.ref)

    // Wait for the JABJEA to be spawned. Then kill it:
    parent.expectNoMessage(max = AwaitAlmostNothing)
    deathwatch.expectNoMessage(max = AwaitAlmostNothing)
    constructionPromise.future onComplete {
      case Success(jabjea) =>
        jabjea ! JabjeaExplode
      case Failure(throwable) =>
        val exception = new RuntimeException("Error creating jabjea for test!", throwable)
        exception.printStackTrace()
        throw exception
    }

    parent.expectMsgPF(max = TimeoutDuration) {
      case JobFailedNonRetryableResponse(_, throwable, _) =>
        throwable.getMessage should be("AwsBatchAsyncBackendJobExecutionActor failed and didn't catch its exception. This condition has been handled and the job will be marked as failed.")
    }
  }
}

class TestAwsBatchJobExecutionActor(params: StandardSyncExecutionActorParams,
                               fakeJabjeaProps: Props) extends StandardSyncExecutionActor(params) {
  override def createAsyncProps(): Props = fakeJabjeaProps
}

class ConstructorFailingJABJEA extends ControllableFailingJabjea(Promise[ActorRef]()) {
  // Explode immediately in the constructor:
  explode()
}

class ControllableFailingJabjea(constructionPromise: Promise[ActorRef]) extends Actor {
  def explode(): Unit = {
    val boom = 1 == 1
    if (boom) throw new RuntimeException("Test Exception! Don't panic if this appears during a test run!")
      with NoStackTrace
  }
  constructionPromise.trySuccess(self)
  override def receive: Receive = {
    case JabjeaExplode => explode()
  }
}

object ControllableFailingJabjea {
  case object JabjeaExplode
}
