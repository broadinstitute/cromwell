package centaur.test

import java.time.OffsetDateTime
import java.util.UUID

import cats.Monad
import centaur._
import centaur.api.CromwellClient._
import centaur.api._
import centaur.test.metadata.WorkflowMetadata
import centaur.test.workflow.Workflow
import scala.annotation.tailrec
import scala.concurrent.{Future, blocking}
import scala.concurrent.duration.FiniteDuration
import scala.util.{Failure, Success, Try}
import scala.concurrent.ExecutionContext.Implicits.global

/**
  * A simplified riff on the final tagless pattern where the interpreter (monad & related bits) are fixed. Operation
  * functions create an instance of a Test and override the run method to do their bidding. It is unlikely that you
  * should be modifying Test directly, instead most likely what you're looking to do is add a function to the Operations
  * object below
  */
sealed abstract class Test[A] {
  def run: Try[A]
}

object Test {
  implicit val testMonad: Monad[Test] = new Monad[Test] {
    override def flatMap[A, B](fa: Test[A])(f: A => Test[B]): Test[B] = {
      new Test[B] {
        override def run: Try[B] = fa.run flatMap { f(_).run }
      }
    }

    override def pure[A](x: A): Test[A] = {
      new Test[A] {
        override def run: Try[A] = Try(x)
      }
    }

    /** Call the default non-stack-safe but correct version of this method. */
    override def tailRecM[A, B](a: A)(f: (A) => Test[Either[A, B]]): Test[B] = defaultTailRecM(a)(f)
  }
}

/**
  * Defines functions which are building blocks for test formulas. Each building block is expected to perform
  * a single task and these tasks can be composed together to form arbitrarily complex test strategies. For instance
  * submitting a workflow, polling until a status is reached, retrieving metadata, verifying some value, delaying for
  * N seconds - these would all be operations.
  *
  * All operations are expected to return a Test type and implement the run method. These can then
  * be composed together via a for comprehension as a test formula and then run by some other entity.
  */
object Operations {
  def submitWorkflow(workflow: Workflow): Test[SubmittedWorkflow] = {
    new Test[SubmittedWorkflow] {
      override def run: Try[SubmittedWorkflow] = CromwellClient.submit(workflow)
    }
  }

  /**
    * Polls until a specific status is reached. If a terminal status which wasn't expected is returned, the polling
    * stops with a failure.
    */
  def pollUntilStatus(workflow: SubmittedWorkflow, expectedStatus: WorkflowStatus): Test[SubmittedWorkflow] = {
    def pollDelay() = blocking { Thread.sleep(10000) } // This could be a lot smarter, including cromwell style backoff
    new Test[SubmittedWorkflow] {
      @tailrec
      def doPerform(allowed404s: Int = 2): SubmittedWorkflow = {
        CromwellClient.status(workflow) match {
          case Success(s) if s == expectedStatus => workflow
          case Success(s: TerminalStatus) => throw new Exception(s"Unexpected terminal status $s but was waiting for $expectedStatus")
          case Failure(f) if f.getMessage.contains("404 Not Found") && allowed404s > 0 =>
            // It's possible that we've started polling prior to the metadata service learning of this workflow
            pollDelay()
            doPerform(allowed404s = allowed404s - 1)
          case Failure(f) => throw f
          case _ =>
            pollDelay()
            doPerform()
        }
      }

      override def run: Try[SubmittedWorkflow] = workflowLengthFutureCompletion(Future { doPerform() })
    }
  }

  def validateMetadata(workflow: SubmittedWorkflow, expectedMetadata: WorkflowMetadata, cacheHitUUID: Option[UUID] = None): Test[Unit] = {
    @tailrec
    def eventually(startTime: OffsetDateTime, timeout: FiniteDuration)(f: => Try[Unit]): Try[Unit] = {
      import scala.concurrent.duration._

      f match {
        case Failure(_) if OffsetDateTime.now().isBefore(startTime.plusSeconds(timeout.toSeconds)) =>
          blocking { Thread.sleep(1.second.toMillis) }
          eventually(startTime, timeout)(f)
        case t => t
      }
    }

    new Test[Unit] {
      def validateMetadata(workflow: SubmittedWorkflow, expectedMetadata: WorkflowMetadata, cacheHitUUID: Option[UUID] = None): Try[Unit] = {
        def checkDiff(diffs: Iterable[String]): Unit = {
          diffs match {
            case d if d.nonEmpty => throw new Exception(s"Invalid metadata response:\n -${d.mkString("\n -")}\n")
            case _ =>
          }
        }

        CromwellClient.metadata(workflow).map(expectedMetadata.diff(_, workflow.id, cacheHitUUID)).map(checkDiff)
      }

      override def run: Try[Unit] = {
        eventually(OffsetDateTime.now(), CentaurConfig.metadataConsistencyTimeout) {
          validateMetadata(workflow, expectedMetadata, cacheHitUUID)
        }
      }
    }
  }

  /**
    * Verify that none of the calls within the workflow are cached.
    */
  def validateNoCacheHits(metadata: WorkflowMetadata, workflowName: String): Test[Unit] = {
    new Test[Unit] {
      override def run: Try[Unit] = {
        val cacheHits = metadata.value.keySet filter { _.contains("Call caching read result") }

        if (cacheHits.isEmpty) Success(())
        else Failure(new Exception(s"Found unexpected cache hits for $workflowName: \n"))
      }
    }
  }

  // FIXME: Should be abstracted w/ validateMetadata - ATM still used by the unused caching tests
  def retrieveMetadata(workflow: SubmittedWorkflow): Test[WorkflowMetadata] = {
    new Test[WorkflowMetadata] {
      override def run: Try[WorkflowMetadata] = CromwellClient.metadata(workflow)
    }
  }

  /* Some enhancements of CromwellApi tools specific to these tests */
  def workflowLengthFutureCompletion[T](x: Future[T]) = awaitFutureCompletion(x, CentaurConfig.maxWorkflowLength)
  def metadataFutureCompletion[T](x: Future[T]) = awaitFutureCompletion(x, CentaurConfig.metadataConsistencyTimeout)
}
