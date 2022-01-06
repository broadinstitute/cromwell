package cromwell.engine.workflow.tokens.large


import akka.actor.{Actor, ActorRef}
import cromwell.core.JobToken.JobTokenType
import cromwell.engine.workflow.tokens.large.LargeScaleJobTokenDispenserActorSpec.RunningJobCounter
import cromwell.engine.workflow.tokens.large.MultipleTokenUsingActor.TokenUsingActorCompletion
import cromwell.engine.workflow.tokens.large.PatientTokenNeedingActor.{AllDone, Begin, ImBusy}

/**
  * Sort of equivalent to a very large workflow -
  *
  * I'm going to spawn a large number of 'job's.
  * I don't expect all of them to be served immediately - so I'll allow them to be patient.
  *
  * Because I'm a good citizen, I'm going to record and return a value representing my "peak concurrent tokens distributed"
  */
class MultipleTokenUsingActor(tokenDispenser: ActorRef, tokenType: JobTokenType, totalJobs: Int, hogGroup: String, globalRunningJobCounter: RunningJobCounter) extends Actor {

  var hasToken: Boolean = false

  var starter: ActorRef = _
  var queueWaits: Seq[Long] = Seq.empty
  var startedJobs = 0
  var runningJobs = 0
  var completedJobs = 0
  var errors: List[String] = List.empty

  var maximumConcurrency = 0

  def receive = {
    case Begin =>
      starter = sender()
      (0 until totalJobs) foreach { i =>
        val jobActor = context.actorOf(PatientTokenNeedingActor.props(tokenDispenser, tokenType, hogGroup), name = self.path.name + s"job$i")
        jobActor ! Begin
      }
      startedJobs = totalJobs
    case ImBusy(msInQueue) =>
      runningJobs += 1
      globalRunningJobCounter.increment
      queueWaits :+= msInQueue
      maximumConcurrency = math.max(runningJobs, maximumConcurrency)
    case AllDone =>
      completedJobs += 1
      runningJobs -= 1
      globalRunningJobCounter.decrement
      if (completedJobs == startedJobs) {
        starter ! TokenUsingActorCompletion(queueWaits, maximumConcurrency, errors)
        context.stop(self)
      }
      ()
    case other =>
      errors :+= s"Received unexpected message $other"
  }
}


object MultipleTokenUsingActor {
  final case class TokenUsingActorCompletion(queueWaits: Seq[Long], maximumConcurrency: Int, errors: List[String])
}
