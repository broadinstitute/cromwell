package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor._
import org.scalatest.{FlatSpec, Matchers}

class EngineJobExecutionActorTransitionsSpec extends FlatSpec with Matchers {

  val callCachingStateCycle = List(
    CheckingCallCache,
    FetchingCachedOutputsFromDatabase,
    BackendIsCopyingCachedOutputs,
    InvalidatingCacheEntry
  )

  "EngineJobExecutionActor transitions" should "not list all cache cycle iterations" in {

    import EngineJobExecutionActorTransitionsSpec.MultipliableList

    val cacheReadCycles = 5

    val longCallCachingCycleStateSequence = List(
      Pending,
      RequestingExecutionToken,
      CheckingJobStore,
      CheckingCallCache,
      FetchingCachedOutputsFromDatabase,
      CheckingCacheEntryExistence) ++ callCachingStateCycle * cacheReadCycles ++ List(
      WaitingForValueStore,
      PreparingJob,
      RunningJob,
      UpdatingCallCache,
      UpdatingJobStore
    )

    longCallCachingCycleStateSequence.length should be(11 + cacheReadCycles * callCachingStateCycle.size)

    val transitionSequence = longCallCachingCycleStateSequence.sliding(2) map {
      case fromState :: toState :: _ => EngineJobExecutionActorState.transitionEventString(fromState, toState)
      case _ => fail("Programmer blunder. This test writer had one job to do...")
    } collect {
      case Some(stateName) => stateName
    }

    transitionSequence.toList should be(List(
      // "Pending", <-- NB: There's no transition into "Pending" because that was the start state
      "RequestingExecutionToken",
      "CheckingJobStore",
      "CallCacheReading",
      "WaitingForValueStore",
      "PreparingJob",
      "RunningJob",
      "UpdatingCallCache",
      "UpdatingJobStore",
    ))
  }
}

object EngineJobExecutionActorTransitionsSpec {
  implicit class MultipliableList[A](val list: List[A]) extends AnyVal {
    final def *(i: Int): List[A] = if (i == 0 ) List.empty else if (i == 1) list else list ++ (list * (i - 1))
  }
}
