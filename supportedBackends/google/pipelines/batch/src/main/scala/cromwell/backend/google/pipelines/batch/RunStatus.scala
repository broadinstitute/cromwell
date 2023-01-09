package cromwell.backend.google.pipelines.batch

//import cromwell.core.ExecutionEvent
//import com.google.cloud.batch.v1.JobStatus
//import scala.util.{Try, Success}

sealed trait RunStatus {
    //import RunStatus._

    //def isTerminal: Boolean

    /*
    def isRunningOrComplete = this match {
        case Running | _: TerminalRunStatus => true
        case _ => false
    }
    */
}

object RunStatus {

    /*
    def fromJobStatus(status: JobStatus.State,  eventList: Seq[ExecutionEvent] = Seq.empty): Try[RunStatus] = {
        status match {
            case JobStatus.State.QUEUED =>
                println("Job is queued")
                Success(Running)
            case JobStatus.State.SCHEDULED =>
                println ("job scheduled")
                Success(Running)
            case JobStatus.State.RUNNING =>
                println("job running")
                Success(Running)
            case JobStatus.State.SUCCEEDED =>
                println("job succeeded")
                //Success(Succeeded(eventList))
                Success(Succeeded)
            case _ =>
                println("no matches in RunStatus")
                Success(Running)
        }

    }
    */

    case object Initializing extends RunStatus {
        //def isTerminal=false
    }
    case object Running extends RunStatus {
        //def isTerminal=false
    }
    case object Succeeded extends RunStatus {
        //def isTerminal=true
    }

    case object TempBatch extends RunStatus {
        //def isTerminal=false
    }
    sealed trait TerminalRunStatus extends RunStatus {
      //def eventList: Seq[ExecutionEvent]
    }

    case class TempBatch()


    //case class Succeeded(eventList: Seq[ExecutionEvent]) extends TerminalRunStatus {

    case class Succeeded() extends TerminalRunStatus {
        //machineType: Option[String],
        //zone: Option[String])
        //instanceName: Option[String]) extends TerminalRunStatus {
        override def toString = "Succeeded"
        def isTerminal = true
    }

}
