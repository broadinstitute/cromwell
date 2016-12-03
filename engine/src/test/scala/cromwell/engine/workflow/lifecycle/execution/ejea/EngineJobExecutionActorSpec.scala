package cromwell.engine.workflow.lifecycle.execution.ejea

import akka.actor.Actor
import akka.testkit.TestFSMRef
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.jobstore.{Pending => _}
import cromwell.CromwellTestKitSpec
import cromwell.backend.BackendJobExecutionActor
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionActorCommand
import cromwell.core.callcaching._
import org.scalatest._
import org.specs2.mock.Mockito

import scala.concurrent.duration._
import scala.language.postfixOps


trait EngineJobExecutionActorSpec extends CromwellTestKitSpec
  with Matchers with Mockito with BeforeAndAfterAll with BeforeAndAfter {

  // If we WANT something to happen, make sure it happens within this window:
  val awaitTimeout: FiniteDuration = 10 seconds
  // If we want something to NOT happen, make sure it doesn't happen within this window. This is almost nothing on
  // the basis that the thing we're awaiting would probably have happened by the time the call was made. Otherwise choose
  // another time!
  val awaitAlmostNothing: FiniteDuration = 100 milliseconds

  implicit override val patienceConfig: PatienceConfig = PatienceConfig(timeout = scaled(awaitTimeout), interval = scaled(awaitAlmostNothing))


  // The default values for these are "null". The helper is created in "before", the ejea is up to the test cases
  private[ejea] var helper: PerTestHelper = _
  private[ejea] var ejea: TestFSMRef[EngineJobExecutionActorState, EJEAData, MockEjea] = _
  implicit def stateUnderTest: EngineJobExecutionActorState

  before {
    helper = new PerTestHelper
  }
  after {
    List(
      ("FetchCachedResultsActor", helper.fetchCachedResultsActorCreations),
      ("JobHashingActor", helper.jobHashingInitializations),
      ("CallCacheInvalidateActor", helper.invalidateCacheActorCreations),
      ("CallCacheWriteActor", helper.callCacheWriteActorCreations)) foreach {
      case (name, GotTooMany(list)) => fail(s"Too many $name creations (${list.size})")
      case _ => // Fine.
    }

    if (ejea != null) system.stop(ejea)
  }

  /**
    *
    * var fetchCachedResultsActorCreations: ExpectOne[(CacheHit, Seq[TaskOutput])] = NothingYet
    * var jobHashingInitializations: ExpectOne[(BackendJobDescriptor, CallCachingActivity)] = NothingYet
    * var callCacheWriteActorCreations: ExpectOne[(CallCacheHashes, SucceededResponse)] = NothingYet
    */

  override def afterAll(): Unit = {
    system.terminate()
    super.afterAll()
  }

  // Some helper lists
  val CallCachingModes = List(CallCachingOff, CallCachingActivity(ReadCache), CallCachingActivity(WriteCache), CallCachingActivity(ReadAndWriteCache))
  case class RestartOrExecuteCommandTuple(operationName: String, restarting: Boolean, expectedMessageToBjea: BackendJobExecutionActorCommand)
  val RestartOrExecuteCommandTuples = List(
    RestartOrExecuteCommandTuple("execute", restarting = false, BackendJobExecutionActor.ExecuteJobCommand),
    RestartOrExecuteCommandTuple("restart", restarting = true, BackendJobExecutionActor.RecoverJobCommand))
}

object EngineJobExecutionActorSpec {
  implicit class EnhancedTestEJEA[S,D,T <: Actor](me: TestFSMRef[S, D, T]) {
    // Like setState, but mirrors back the EJEA (for easier inlining)
    def setStateInline(state: S = me.stateName, data: D = me.stateData) = {
      me.setState(state, data)
      me
    }
  }
}


