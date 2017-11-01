package centaur.test

import java.time.OffsetDateTime
import java.util.UUID

import cats.Monad
import centaur._
import centaur.api.CentaurCromwellClient
import centaur.api.CentaurCromwellClient.sendReceiveFutureCompletion
import centaur.test.metadata.WorkflowMetadata
import centaur.test.submit.SubmitResponse
import centaur.test.workflow.Workflow
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model.Operation
import com.google.auth.http.HttpCredentialsAdapter
import com.google.auth.oauth2.GoogleCredentials
import com.google.cloud.compute.{ComputeOptions, InstanceId}
import cromwell.api.CromwellClient.UnsuccessfulRequestException
import cromwell.api.model.{Failed, SubmittedWorkflow, TerminalStatus, WorkflowStatus}
import spray.json.JsString

import scala.annotation.tailrec
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration.FiniteDuration
import scala.concurrent.{Await, Future, blocking}
import scala.util.{Failure, Success, Try}

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
  def successful[A](value: A): Test[A] = testMonad.pure(value)
  def failed[A](exception: Exception) = new Test[A] {
    override def run = Failure(exception)
  }

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
    override def tailRecM[A, B](a: A)(f: (A) => Test[Either[A, B]]): Test[B] = {
      flatMap(f(a)) {
        case Right(b) => pure(b)
        case Left(nextA) => tailRecM(nextA)(f)
      }
    }
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
  lazy val jsonFactory = JacksonFactory.getDefaultInstance
  lazy val httpTransport = GoogleNetHttpTransport.newTrustedTransport
  lazy val genomics = new Genomics.Builder(
    httpTransport,
    jsonFactory,
    new HttpCredentialsAdapter(GoogleCredentials.getApplicationDefault)
  ).setApplicationName("centaur")
    .setRootUrl("https://genomics.googleapis.com/")
    .build()

  val compute = ComputeOptions.getDefaultInstance.getService

  def submitWorkflow(workflow: Workflow): Test[SubmittedWorkflow] = {
    new Test[SubmittedWorkflow] {
      override def run: Try[SubmittedWorkflow] = CentaurCromwellClient.submit(workflow)
    }
  }

  def submitInvalidWorkflow(workflow: Workflow): Test[SubmitResponse] = {
    new Test[SubmitResponse] {
      override def run: Try[SubmitResponse] = {
        Try {
          CentaurCromwellClient.submit(workflow) match {

            case Success(submittedWorkflow) =>
              throw new RuntimeException(
                s"Expected a failure but got a successfully submitted workflow with id ${submittedWorkflow.id}")

            case Failure(unsuccessfulRequestException: UnsuccessfulRequestException) =>
              val httpResponse = unsuccessfulRequestException.httpResponse
              val statusCode = httpResponse.status.intValue()
              val message = httpResponse.entity match {
                case akka.http.scaladsl.model.HttpEntity.Strict(_, data) => data.utf8String
                case _ =>
                  throw new RuntimeException(s"Expected a strict http response entity but got ${httpResponse.entity}")
              }
              SubmitResponse(statusCode, message)

            case Failure(unexpected) => throw unexpected
          }
        }
      }
    }
  }

  def abortWorkflow(workflow: SubmittedWorkflow) = {
    new Test[WorkflowStatus] {
      override def run: Try[WorkflowStatus] = CentaurCromwellClient.abort(workflow)
    }
  }
  
  def waitFor(duration: FiniteDuration) = {
    import scala.concurrent.duration._

    new Test[Unit] {
      override def run = {
        // Give some margin to the timeout
        Try(Await.result(Future { Thread.sleep(duration.toMillis) }, duration.plus(2.seconds)))
      }
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
        CentaurCromwellClient.status(workflow) match {
          case Success(s) if s == expectedStatus => workflow
          case Success(s: TerminalStatus) => throw new Exception(s"Unexpected terminal status $s but was waiting for $expectedStatus")
          case Failure(f) if f.getMessage.contains("404 Not Found") && allowed404s > 0 =>
            // It's possible that we've started polling prior to the metadata service learning of this workflow
            pollDelay()
            doPerform(allowed404s = allowed404s - 1)
          case Failure(f) if CromwellManager.isReady => throw f
          case _ =>
            pollDelay()
            doPerform()
        }
      }

      override def run: Try[SubmittedWorkflow] = workflowLengthFutureCompletion(() => Future { doPerform() })
    }
  }

  /**
    * Validate that the given jobId matches the one in the metadata
    */
  def validateRecovered(workflow: SubmittedWorkflow, callFqn: String, formerJobId: String): Test[Unit] = {
    new Test[Unit] {
      def doPerform(): Unit = {
        CentaurCromwellClient.metadata(workflow) match {
          case Success(s) =>
            s.value.get(s"calls.$callFqn.jobId") match {
              case Some(newJobId) if newJobId.asInstanceOf[JsString].value == formerJobId => ()
              case Some(_) => throw new Exception("Pre-restart job ID did not match post restart job ID")
              case _ => throw new Exception("Cannot find a post restart job ID")
            }
          case Failure(f) => throw f
        }
      }

      override def run: Try[Unit] = sendReceiveFutureCompletion(() => Future { doPerform() })
    }
  }

  def validatePAPIAborted(jobId: String, workflow: SubmittedWorkflow): Test[Unit] = {
    import scala.concurrent.duration._

    def eventually(startTime: OffsetDateTime, timeout: FiniteDuration)(f: => Try[Unit]): Try[Unit] = {
      f match {
        case Failure(_) if OffsetDateTime.now().isBefore(startTime.plusSeconds(timeout.toSeconds)) =>
          blocking { Thread.sleep(1.second.toMillis) }
          eventually(startTime, timeout)(f)
        case t => t
      }
    }
    
    new Test[Unit] {
      def checkVMTerminated(): Unit = {
          CentaurCromwellClient.metadata(workflow) match {
            case Success(metadata) =>
              // Note: this doesn't work as of now because Cromwell doesn't publish metadata values for instanceName and
              // zone if the job gets cancelled
              val instanceName = metadata.value.collectFirst({
                case (key, value) if key.endsWith("jes.instanceName") => value.asInstanceOf[JsString].value
              }).getOrElse(throw new Exception("Cannot find the instance name in metadata"))
              
              val zone = metadata.value.collectFirst({
                case (key, value) if key.endsWith("jes.zone") => value.asInstanceOf[JsString].value
              }).getOrElse(throw new Exception("Cannot find the zone in metadata"))
              
              val instanceId = InstanceId.of(zone, instanceName)

              Try(compute.getInstance(instanceId)) match {
                case Success(null) => // all good, couldn't find the instance
                case Success(_) => throw new Exception("Underlying VM is still running")
                case Failure(f) => throw f
              }
            case Failure(f) => throw f
          }
      }
      
      def checkPAPIAborted(): Unit = {
        val operation: Operation = genomics.operations().get(jobId).execute()
        val done = operation.getDone
        val aborted = operation.getError.getCode == 1 && operation.getError.getMessage.startsWith("Operation canceled")
        if (!(done && aborted)) {
          throw new Exception("Underlying JES job was not aborted properly")
        }
      }

      override def run: Try[Unit] = if (jobId.startsWith("operations/")) {
        // The PAPI status should be aborted immediately
        // Note: this doesn't work as of now because Cromwell considers the workflow aborted
        // as soon as it has requested cancellation, not when PAPI says its cancelled
        Try(checkPAPIAborted())
        // Give some time to the VM to actually die (PAPI says it's cancelled before the VM is actually killed)
        eventually(OffsetDateTime.now(), 1.minute) {
          Try(checkVMTerminated())
        }
      } else Success(())
    }
  }

  /**
    * Polls until a specific call is in Running state. Returns the job id.
    */
  def pollUntilCallIsRunning(workflow: SubmittedWorkflow, callFqn: String): Test[String] = {
    // We want to keep this smaller than the runtime of the call we're polling for
    // For JES it should be fine but locally it can be quite fast
    def pollDelay() = blocking { Thread.sleep(5000) }

    def findCallStatus(metadata: WorkflowMetadata): Option[(String, String)] = {
      for {
        status <- metadata.value.get(s"calls.$callFqn.executionStatus")
        jobId <- metadata.value.get(s"calls.$callFqn.jobId")
      } yield (status.asInstanceOf[JsString].value, jobId.asInstanceOf[JsString].value)
    }

    new Test[String] {
      @tailrec
      def doPerform(allowed404s: Int = 2): String = {
        val metadata = for {
        // We don't want to keep going forever if the workflow failed
          status <- CentaurCromwellClient.status(workflow)
          _ <- status match {
            case Failed => Failure(new Exception("Workflow Failed"))
            case _ => Success(())
          }
          metadata <- CentaurCromwellClient.metadata(workflow)
        } yield metadata

        metadata match {
          case Success(s) =>
            findCallStatus(s) match {
              case Some(("Running", jobId)) => jobId
              case Some(("Failed", _)) => throw new Exception(s"$callFqn failed")
              case _ =>
                pollDelay()
                doPerform()
            }
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

      override def run: Try[String] = workflowLengthFutureCompletion(() => Future { doPerform() })
    }
  }

  def validateMetadata(submittedWorkflow: SubmittedWorkflow, workflowSpec: Workflow, cacheHitUUID: Option[UUID] = None): Test[WorkflowMetadata] = {
    @tailrec
    def eventually(startTime: OffsetDateTime, timeout: FiniteDuration)(f: => Try[WorkflowMetadata]): Try[WorkflowMetadata] = {
      import scala.concurrent.duration._

      f match {
        case Failure(_) if OffsetDateTime.now().isBefore(startTime.plusSeconds(timeout.toSeconds)) =>
          blocking { Thread.sleep(1.second.toMillis) }
          eventually(startTime, timeout)(f)
        case t => t
      }
    }

    new Test[WorkflowMetadata] {
      def validateMetadata(workflow: SubmittedWorkflow, expectedMetadata: WorkflowMetadata, cacheHitUUID: Option[UUID] = None): Try[WorkflowMetadata] = {
        def checkDiff(diffs: Iterable[String]): Unit = {
          diffs match {
            case d if d.nonEmpty => throw new Exception(s"Invalid metadata response:\n -${d.mkString("\n -")}\n")
            case _ =>
          }
        }
        cleanUpImports(workflow)
        
        def validateUnwantedMetadata(actualMetadata: WorkflowMetadata) = if (workflowSpec.notInMetadata.nonEmpty) {
          // Check that none of the "notInMetadata" keys are in the actual metadata
          val absentMdDif = workflowSpec.notInMetadata.toSet.diff(actualMetadata.value.keySet)
          if (absentMdDif.isEmpty) throw new Exception(s"Found unwanted keys in metadata: ${absentMdDif.mkString(", ")}")
        }

        for {
          actualMetadata <- CentaurCromwellClient.metadata(workflow)
          _ = validateUnwantedMetadata(actualMetadata)
          diffs = expectedMetadata.diff(actualMetadata, workflow.id.id, cacheHitUUID)
          _ = checkDiff(diffs)
        } yield actualMetadata
      }

      override def run: Try[WorkflowMetadata] = workflowSpec.metadata match {
        case Some(expectedMetadata) =>
          eventually(OffsetDateTime.now(), CentaurConfig.metadataConsistencyTimeout) {
            validateMetadata(submittedWorkflow, expectedMetadata, cacheHitUUID)
          }
        // Nothing to wait for, so just return the first metadata we get back:
        case None => CentaurCromwellClient.metadata(submittedWorkflow)
      }
    }
  }

  /**
    * Verify that none of the calls within the workflow are cached.
    */
  def validateCacheResultField(metadata: WorkflowMetadata, workflowName: String, blacklistedValue: String): Test[Unit] = {
    new Test[Unit] {
      override def run: Try[Unit] = {
        val badCacheResults = metadata.value collect {
          case (k, JsString(v)) if k.contains("callCaching.result") && v.contains(blacklistedValue) => s"$k: $v"
        }

        if (badCacheResults.isEmpty) Success(())
        else Failure(new Exception(s"Found unexpected cache hits for $workflowName:${badCacheResults.mkString("\n", "\n", "\n")}"))
      }
    }
  }

  def validateDirectoryContentsCounts(workflowDefinition: Workflow, submittedWorkflow: SubmittedWorkflow): Test[Unit] = new Test[Unit] {
    private val workflowId = submittedWorkflow.id.id.toString
    override def run: Try[Unit] = workflowDefinition.directoryContentCounts match {
      case None => Success(())
      case Some(directoryContentCountCheck) =>
        val counts = directoryContentCountCheck.expectedDrectoryContentsCounts map {
          case (directory, count) =>
            val substitutedDir = directory.replaceAll("<<UUID>>", workflowId)
            (substitutedDir, count, directoryContentCountCheck.checkFiles.countObjectsAtPath(substitutedDir))
        }

        val badCounts = counts collect {
          case (directory, expectedCount, actualCount) if expectedCount != actualCount => s"Expected to find $expectedCount item(s) at $directory but got $actualCount"
        }
        if (badCounts.isEmpty) Success(()) else Failure(new Exception(badCounts.mkString("\n", "\n", "\n")))
    }
  }

  def validateNoCacheHits(metadata: WorkflowMetadata, workflowName: String): Test[Unit] = validateCacheResultField(metadata, workflowName, "Cache Hit")
  def validateNoCacheMisses(metadata: WorkflowMetadata, workflowName: String): Test[Unit] = validateCacheResultField(metadata, workflowName, "Cache Miss")

  def validateSubmitFailure(workflow: Workflow,
                            expectedSubmitResponse: SubmitResponse,
                            actualSubmitResponse: SubmitResponse): Test[Unit] = {
    new Test[Unit] {
      override def run: Try[Unit] = {
        Try {
          if (expectedSubmitResponse == actualSubmitResponse) {
            ()
          } else {
            throw new RuntimeException(
              s"""|
                  |Expected
                  |$expectedSubmitResponse
                  |
                  |but got:
                  |$actualSubmitResponse
                  |""".stripMargin
            )
          }
        }
      }
    }
  }

  /**
    * Clean up temporary zip files created for Imports testing.
    */
  def cleanUpImports(submittedWF: SubmittedWorkflow) = {
    submittedWF.workflow.zippedImports match {
      case Some(zipFile) => zipFile.delete(swallowIOExceptions = true)
      case None => //
    }
  }

  // FIXME: Should be abstracted w/ validateMetadata - ATM still used by the unused caching tests
  def retrieveMetadata(workflow: SubmittedWorkflow): Test[WorkflowMetadata] = {
    new Test[WorkflowMetadata] {
      override def run: Try[WorkflowMetadata] = CentaurCromwellClient.metadata(workflow)
    }
  }

  /* Some enhancements of CromwellApi tools specific to these tests */
  def workflowLengthFutureCompletion[T](x: () => Future[T]) = CentaurCromwellClient.maxWorkflowLengthCompletion(x)
}
