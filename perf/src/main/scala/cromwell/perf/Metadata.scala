package cromwell.perf

import java.time.{Duration, OffsetDateTime}

case class Metadata(id: String,
                    workflowName: String,
                    submission: OffsetDateTime,
                    start: OffsetDateTime,
                    end: OffsetDateTime,
                    status: String,
                    calls: Option[Map[String, Seq[Call]]]) {

  private def addInt(x: Int, y: Int): Int = x + y

  private def addDuration(x: Duration, y: Duration): Duration = x.plus(y)

  private def sumElementsInOptionSeq[A](listOption: Option[Iterable[A]], op: (A, A) => A, default: A): A = {
    listOption match {
      case Some(list) => list.reduce[A]((a, b) => op(a,b))
      case None => default
    }
  }

  /**
    * @return Time between the submission and start of workflow
    */
  val workflowStartedAfter: Duration = Duration.between(submission, start)

  val workflowRunningTime: Duration = Duration.between(start, end)

  val totalJobsPerRootWf: Int = sumElementsInOptionSeq(calls.map(taskMap => taskMap.map(callsPerTask => callsPerTask._2.size)), addInt, 0)

  val avgCacheRetries: Int = {
    if (totalJobsPerRootWf > 0) {
      val cacheRetriesList = calls.map(taskMap => taskMap.flatMap(callsPerTask => callsPerTask._2.map(call => call.cacheCopyRetries)))
      sumElementsInOptionSeq(cacheRetriesList, addInt, 0) / totalJobsPerRootWf
    }
    else 0
  }

  val avgTimeInCallCachingState: Duration = {
    if (totalJobsPerRootWf > 0) {
      val timeInCallCachingStateList = calls.map(taskMap => taskMap.flatMap(callsPerTask => callsPerTask._2.map(call => call.timeInCallCachingState)))
      sumElementsInOptionSeq(timeInCallCachingStateList, addDuration, Duration.ZERO).dividedBy(totalJobsPerRootWf.toLong)
    }
    else Duration.ZERO
  }

  val avgTimeInCheckingCallCacheState: Duration = {
    if (totalJobsPerRootWf > 0) {
      val timeInCheckingCallCacheStateList = calls.map(taskMap => taskMap.flatMap(callsPerTask => callsPerTask._2.map(call => call.timeInCheckingCallCacheState)))
      sumElementsInOptionSeq(timeInCheckingCallCacheStateList, addDuration, Duration.ZERO).dividedBy(totalJobsPerRootWf.toLong)
    }
    else Duration.ZERO
  }

  val avgTimeInJobPreparation: Duration = {
    if (totalJobsPerRootWf > 0) {
      val timeInJobPreparationList = calls.map(taskMap => taskMap.flatMap(callsPerTask => callsPerTask._2.map(call => call.timeInJobPreparation)))
      sumElementsInOptionSeq(timeInJobPreparationList, addDuration, Duration.ZERO).dividedBy(totalJobsPerRootWf.toLong)
    }
    else Duration.ZERO
  }

  val avgTimeForFetchingAndCopyingCacheHit: Duration = {
    if (totalJobsPerRootWf > 0) {
      val timeInCopyingList = calls.map(taskMap => taskMap.flatMap(callsPerTask => callsPerTask._2.map(call => call.timeForFetchingAndCopyingCacheHit)))
      sumElementsInOptionSeq(timeInCopyingList, addDuration, Duration.ZERO).dividedBy(totalJobsPerRootWf.toLong)
    }
    else Duration.ZERO
  }
}


