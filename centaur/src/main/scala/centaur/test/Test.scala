package centaur.test

import java.util.UUID

import cats.Monad
import cats.effect.IO
import cats.instances.list._
import cats.syntax.traverse._
import centaur._
import centaur.api.CentaurCromwellClient
import centaur.test.metadata.WorkflowFlatMetadata
import centaur.test.metadata.WorkflowFlatMetadata._
import centaur.test.submit.SubmitHttpResponse
import centaur.test.workflow.Workflow
import com.google.api.services.genomics.{Genomics, GenomicsScopes}
import com.google.api.services.storage.StorageScopes
import com.google.auth.Credentials
import com.google.auth.http.HttpCredentialsAdapter
import com.google.auth.oauth2.ServiceAccountCredentials
import com.google.cloud.storage.{Storage, StorageOptions}
import com.typesafe.config.Config
import com.typesafe.scalalogging.StrictLogging
import common.validation.Validation._
import configs.syntax._
import cromwell.api.CromwellClient.UnsuccessfulRequestException
import cromwell.api.model.{CallCacheDiff, Failed, SubmittedWorkflow, Succeeded, TerminalStatus, WaasDescription, WorkflowId, WorkflowMetadata, WorkflowStatus}
import cromwell.cloudsupport.aws.AwsConfiguration
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import io.circe.parser._
import mouse.all._
import org.apache.commons.lang3.exception.ExceptionUtils
import software.amazon.awssdk.auth.credentials.{AwsBasicCredentials, StaticCredentialsProvider}
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.s3.S3Client
import spray.json._

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.TimeoutException
import scala.concurrent.duration._
import scala.util.{Failure, Try}

/**
  * A simplified riff on the final tagless pattern where the interpreter (monad & related bits) are fixed. Operation
  * functions create an instance of a Test and override the run method to do their bidding. It is unlikely that you
  * should be modifying Test directly, instead most likely what you're looking to do is add a function to the Operations
  * object below
  */
sealed abstract class Test[A] {
  def run: IO[A]
}

object Test {
  def successful[A](value: A): Test[A] = testMonad.pure(value)

  def invalidTestDefinition[A](message: String, workflowDefinition: Workflow): Test[A] = {
    new Test[A] {
      override def run: IO[Nothing] = IO.raiseError(CentaurTestException(message, workflowDefinition))
    }
  }

  implicit val testMonad: Monad[Test] = new Monad[Test] {
    override def flatMap[A, B](fa: Test[A])(f: A => Test[B]): Test[B] = {
      new Test[B] {
        override def run: IO[B] = fa.run flatMap { f(_).run }
      }
    }

    override def pure[A](x: A): Test[A] = {
      new Test[A] {
        override def run: IO[A] = IO.pure(x)
      }
    }

    /** Call the default non-stack-safe but correct version of this method. */
    override def tailRecM[A, B](a: A)(f: A => Test[Either[A, B]]): Test[B] = {
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
object Operations extends StrictLogging {
  lazy val configuration: GoogleConfiguration = GoogleConfiguration(CentaurConfig.conf)
  lazy val googleConf: Config = CentaurConfig.conf.getConfig("google")
  lazy val authName: String = googleConf.getString("auth")
  lazy val genomicsEndpointUrl: String = googleConf.getString("genomics.endpoint-url")
  lazy val genomicsAndStorageScopes = List(StorageScopes.CLOUD_PLATFORM_READ_ONLY, GenomicsScopes.GENOMICS)
  lazy val credentials: Credentials = configuration.auth(authName)
    .unsafe
    .credentials(genomicsAndStorageScopes)
  lazy val credentialsProjectOption: Option[String] = {
    Option(credentials) collect {
      case serviceAccountCredentials: ServiceAccountCredentials => serviceAccountCredentials.getProjectId
    }
  }
  lazy val confProjectOption: Option[String] = googleConf.get[Option[String]]("project") valueOrElse None
  // The project from the config or from the credentials. By default the project is read from the system environment.
  lazy val projectOption: Option[String] = confProjectOption orElse credentialsProjectOption

  lazy val genomics: Genomics = {
    val builder = new Genomics.Builder(
      GoogleAuthMode.httpTransport,
      GoogleAuthMode.jsonFactory,
      new HttpCredentialsAdapter(credentials)
    )
    builder
      .setApplicationName(configuration.applicationName)
      .setRootUrl(genomicsEndpointUrl)
      .build()
  }

  lazy val storage: Storage = {
    val builder = StorageOptions.newBuilder().setCredentials(credentials)
    projectOption foreach builder.setProjectId
    val storageOptions = builder.build()
    storageOptions.getService
  }

  implicit private val timer = IO.timer(global)
  implicit private val contextShift = IO.contextShift(global)

  lazy val awsConfiguration: AwsConfiguration = AwsConfiguration(CentaurConfig.conf)
  lazy val awsConf: Config = CentaurConfig.conf.getConfig("aws")
  lazy val awsAuthName: String = awsConf.getString("auths")
  lazy val region: String  = awsConf.getString("region")
  lazy val accessKeyId: String  = awsConf.getString("access-key")
  lazy val secretAccessKey: String = awsConf.getString("secret-key")

  def buildAmazonS3Client: S3Client = {
    val basicAWSCredentials = AwsBasicCredentials.create(accessKeyId, secretAccessKey)
    S3Client.builder()
      .region(Region.of(region))
      .credentialsProvider(StaticCredentialsProvider.create(basicAWSCredentials))
      .build()
  }

  def submitWorkflow(workflow: Workflow): Test[SubmittedWorkflow] = {
    new Test[SubmittedWorkflow] {
      override def run: IO[SubmittedWorkflow] = for {
        id <- CentaurCromwellClient.submit(workflow)
      } yield id
    }
  }

  /**
    * A smoke test of the version endpoint, this confirms that a) nothing explodes and b) the result must be a JSON object
    * with a "cromwell" key. There is no sanity checking of the version value here.
    */
  def checkVersion(): Test[Unit] = new Test[Unit] {
    override def run: IO[Unit] = for {
      _ <- CentaurCromwellClient.version
    } yield ()
  }

  def checkDescription(workflow: Workflow, validityExpectation: Option[Boolean], retries: Int = 3): Test[Unit] = {
    new Test[Unit] {

      val timeout = 60.seconds

      def checkDescriptionInner(alreadyTried: Int): IO[Unit] = {
        val timeoutStackTraceString = ExceptionUtils.getStackTrace(new Exception)
        (CentaurCromwellClient.describe(workflow) flatMap { d: WaasDescription =>
          validityExpectation match {
            case None => IO.pure(())
            case Some(d.valid) => IO.pure(())
            case Some(otherExpectation) =>
              logger.error(s"Unexpected 'valid=${d.valid}' response when expecting $otherExpectation. Full unexpected description:${System.lineSeparator()}$d")
              IO.raiseError(new Exception(s"Expected this workflow's /describe validity to be '$otherExpectation' but got: '${d.valid}' (errors: ${d.errors.mkString(", ")})"))
          }
        }).timeoutTo(timeout, IO {
          if (alreadyTried + 1 >= retries) {
            throw new TimeoutException("Timeout from checkDescription 60 seconds: " + timeoutStackTraceString)
          } else {
            logger.warn(s"checkDescription timeout on attempt ${alreadyTried + 1}. ")
            checkDescriptionInner(alreadyTried + 1)
            ()
          }
        })
      }


      override def run: IO[Unit] = {


        // We can't describe workflows based on zipped imports, so don't try:
        if (workflow.skipDescribeEndpointValidation || workflow.data.zippedImports.nonEmpty) {
          IO.pure(())
        } else {
          checkDescriptionInner(0)
        }
      }
    }
  }

  def submitInvalidWorkflow(workflow: Workflow): Test[SubmitHttpResponse] = {
    new Test[SubmitHttpResponse] {
      override def run: IO[SubmitHttpResponse] = {
        CentaurCromwellClient.submit(workflow).redeemWith(
          {
            case unsuccessfulRequestException: UnsuccessfulRequestException =>
              val httpResponse = unsuccessfulRequestException.httpResponse
              val statusCode = httpResponse.status.intValue()
              httpResponse.entity match {
                case akka.http.scaladsl.model.HttpEntity.Strict(_, data) =>
                  IO.pure(SubmitHttpResponse(statusCode, data.utf8String))
                case _ =>
                  val message = s"Expected a strict http response entity but got ${httpResponse.entity}"
                  IO.raiseError(CentaurTestException(message, workflow, unsuccessfulRequestException))
              }
            case unexpected: Exception =>
              val message = s"Unexpected error: ${unexpected.getMessage}"
              IO.raiseError(CentaurTestException(message, workflow, unexpected))
            case throwable: Throwable => throw throwable
          },
          {
            submittedWorkflow => {
              val message =
                s"Expected a failure but got a successfully submitted workflow with id ${submittedWorkflow.id}"
              IO.raiseError(CentaurTestException(message, workflow))
            }
          }
        )
      }
    }
  }

  def abortWorkflow(workflow: SubmittedWorkflow) = {
    new Test[WorkflowStatus] {
      override def run: IO[WorkflowStatus] = CentaurCromwellClient.abort(workflow)
    }
  }

  def waitFor(duration: FiniteDuration) = {
    new Test[Unit] {
      override def run = IO.sleep(duration)
    }
  }

  /**
    * Polls until a valid status is reached.
    * If an unexpected terminal status is returned, the polling stops with a failure.
    * If none of the expected statuses are returned in time, the polling also stops with a failure.
    */
  def expectSomeProgress(workflow: SubmittedWorkflow,
                         testDefinition: Workflow,
                         expectedStatuses: Set[WorkflowStatus],
                         timeout: FiniteDuration): Test[SubmittedWorkflow] = {
    new Test[SubmittedWorkflow] {
      def status(remainingTimeout: FiniteDuration): IO[SubmittedWorkflow] = {
        for {
          workflowStatus <- CentaurCromwellClient.status(workflow)
          mappedStatus <- workflowStatus match {
            case s if expectedStatuses.contains(s) => IO.pure(workflow)
            case s: TerminalStatus =>
              CentaurCromwellClient.metadata(workflow) flatMap { metadata =>
                val message = s"Unexpected terminal status $s while waiting for one of [${expectedStatuses.mkString(", ")}] (workflow ID: ${workflow.id})"
                IO.raiseError(CentaurTestException(message, testDefinition, workflow, metadata))
              }
            case _ if remainingTimeout > 0.seconds =>
              for {
                _ <- IO.sleep(10.seconds)
                s <- status(remainingTimeout - 10.seconds)
              } yield s
            case other =>
              val message = s"Cromwell failed to progress into any of the statuses [${expectedStatuses.mkString(", ")}]. Was still '$other' after $timeout (workflow ID: ${workflow.id})"
              IO.raiseError(CentaurTestException(message, testDefinition, workflow))
          }
        } yield mappedStatus
      }

      override def run: IO[SubmittedWorkflow] = status(timeout).timeout(CentaurConfig.maxWorkflowLength)

    }
  }

  /**
    * Polls until a specific status is reached.
    * If a terminal status which wasn't expected is returned, the polling stops with a failure.
    * If the workflow does not become either (a) Running or (b) terminal within 1 minute, the polling stops with a failure.
    */
  def pollUntilStatus(workflow: SubmittedWorkflow,
                      testDefinition: Workflow,
                      expectedStatus: WorkflowStatus): Test[SubmittedWorkflow] = {
    new Test[SubmittedWorkflow] {
      def status: IO[SubmittedWorkflow] = {
        for {
          workflowStatus <- CentaurCromwellClient.status(workflow)
          mappedStatus <- workflowStatus match {
            case s if s == expectedStatus => IO.pure(workflow)
            case s: TerminalStatus =>
              CentaurCromwellClient.metadata(workflow) flatMap { metadata =>
                val failuresString = if (expectedStatus == Succeeded) {
                  (for {
                    metadataJson <- parse(metadata.value).toOption
                    asObject <- metadataJson.asObject
                    failures <- asObject.toMap.get("failures")
                  } yield s" Metadata 'failures' content: ${failures.spaces2}").getOrElse("No additional failure information found in metadata.")
                } else {
                  ""
                }

                val message = s"Unexpected terminal status $s but was waiting for $expectedStatus (workflow ID: ${workflow.id}).$failuresString"
                IO.raiseError(CentaurTestException(message, testDefinition, workflow, metadata))
              }
            case _ => for {
              _ <- IO.sleep(10.seconds)
              s <- status
            } yield s
          }
        } yield mappedStatus
      }

      override def run: IO[SubmittedWorkflow] = status.timeout(CentaurConfig.maxWorkflowLength)
    }
  }

  /**
    * Validate that the given jobId matches the one in the metadata
    */
  def validateRecovered(workflowDefinition: Workflow,
                        workflow: SubmittedWorkflow,
                        metadata: WorkflowMetadata,
                        callFqn: String,
                        formerJobId: String): Test[Unit] = {
    new Test[Unit] {
      override def run: IO[Unit] = CentaurCromwellClient.metadata(workflow) flatMap { s =>
        s.asFlat.value.get(s"calls.$callFqn.jobId") match {
          case Some(newJobId) if newJobId.asInstanceOf[JsString].value == formerJobId => IO.unit
          case Some(newJobId) =>
            val message = s"Pre-restart job ID $formerJobId did not match post restart job ID $newJobId"
            IO.raiseError(CentaurTestException(message, workflowDefinition, workflow, metadata))
          case _ =>
            val message = s"Cannot find a post restart job ID to match pre-restart job ID $formerJobId"
            IO.raiseError(CentaurTestException(message, workflowDefinition, workflow, metadata))
        }
      }
    }
  }

  def validatePAPIAborted(workflowDefinition: Workflow, workflow: SubmittedWorkflow, jobId: String): Test[Unit] = {
    new Test[Unit] {
      def checkPAPIAborted(): IO[Unit] = {
        for {
          operation <- IO { genomics.operations().get(jobId).execute() }
          done = operation.getDone
          operationError = Option(operation.getError)
          aborted = operationError.exists(_.getCode == 1) && operationError.exists(_.getMessage.startsWith("Operation canceled"))
          result <- if (!(done && aborted)) {
            CentaurCromwellClient.metadata(workflow) flatMap { metadata =>
              val message = s"Underlying JES job was not aborted properly. " +
                s"Done = $done. Error = ${operationError.map(_.getMessage).getOrElse("N/A")} (workflow ID: ${workflow.id})"
              IO.raiseError(CentaurTestException(message, workflowDefinition, workflow, metadata))
            }
          } else IO.unit
        } yield result
      }

      override def run: IO[Unit] = if (jobId.startsWith("operations/")) {
        checkPAPIAborted()
      } else IO.unit
    }
  }

  /**
    * Polls until a specific call is in Running state. Returns the job id.
    */
  def pollUntilCallIsRunning(workflowDefinition: Workflow, workflow: SubmittedWorkflow, callFqn: String): Test[String] = {
    // Special case for sub workflow testing
    def findJobIdInSubWorkflow(subWorkflowId: String): IO[Option[String]] = {
      for {
        metadata <- CentaurCromwellClient
          .metadataWithId(WorkflowId.fromString(subWorkflowId))
          .redeem(_ => None, Option.apply)
        jobId <- IO.pure(metadata.flatMap(_.asFlat.value.get("calls.inner_abort.aborted.jobId")))
      } yield jobId.map(_.asInstanceOf[JsString].value)
    }

    def valueAsString(key: String, metadata: WorkflowMetadata) = {
      metadata.asFlat.value.get(key).map(_.asInstanceOf[JsString].value)
    }

    def findCallStatus(metadata: WorkflowMetadata): IO[Option[(String, String)]] = {
      val status = metadata.asFlat.value.get(s"calls.$callFqn.executionStatus")
      val statusString = status.map(_.asInstanceOf[JsString].value)

      for {
        jobId <- valueAsString(s"calls.$callFqn.jobId", metadata).map(jobId => IO.pure(Option(jobId)))
          .orElse(valueAsString(s"calls.$callFqn.subWorkflowId", metadata).map(findJobIdInSubWorkflow))
          .getOrElse(IO.pure(None))
        pair = (statusString, jobId) match {
          case (Some(s), Some(j)) => Option(s -> j)
          case _ => None
        }
      } yield pair
    }

    new Test[String] {
      def doPerform(): IO[String] = {
        for {
          // We don't want to keep going forever if the workflow failed
          status <- CentaurCromwellClient.status(workflow)
          metadata <- CentaurCromwellClient.metadata(workflow)
          _ <- status match {
            case Failed =>
              val message = s"$callFqn failed"
              IO.raiseError(CentaurTestException(message, workflowDefinition, workflow, metadata))
            case _ => IO.unit
          }
          callStatus <- findCallStatus(metadata)
          result <- callStatus match {
            case Some(("Running", jobId)) => IO.pure(jobId)
            case Some(("Failed", _)) =>
              val message = s"$callFqn failed"
              IO.raiseError(CentaurTestException(message, workflowDefinition, workflow, metadata))
            case _ => for {
              _ <- IO.sleep(5.seconds)
              recurse <- doPerform()
            } yield recurse
          }
        } yield result
      }

      override def run: IO[String] = doPerform().timeout(CentaurConfig.maxWorkflowLength)
    }
  }

  def printHashDifferential(workflowA: SubmittedWorkflow, workflowB: SubmittedWorkflow) = new Test[Unit] {
    def hashDiffOfAllCalls = {
      // Extract the workflow name followed by call name to use in the call cache diff endpoint
      val callNameRegexp = """calls\.([^.]*\.[^.]*)\..*""".r

      for {
        md <- CentaurCromwellClient.metadata(workflowB)
        calls = md.asFlat.value.keySet.flatMap({
          case callNameRegexp(name) => Option(name)
          case _ => None
        })
        diffs <- calls.toList.traverse[IO, CallCacheDiff]({ callName =>
          CentaurCromwellClient.callCacheDiff(workflowA, callName, workflowB, callName)
        })
      } yield diffs.flatMap(_.hashDifferential)
    }

    override def run = {
      hashDiffOfAllCalls map {
        case diffs if diffs.nonEmpty && CentaurCromwellClient.LogFailures =>
          Console.err.println(s"Hash differential for ${workflowA.id} and ${workflowB.id}")
          diffs.map({ diff =>
            s"For key ${diff.hashKey}:\nCall A: ${diff.callA.getOrElse("N/A")}\nCall B: ${diff.callB.getOrElse("N/A")}"
          }).foreach(Console.err.println)
        case _ =>
      }
    }
  }

  /* Select only those flat metadata items whose keys begin with the specified prefix, removing the prefix from the keys. Also
   * perform variable substitutions for UUID and WORKFLOW_ROOT and remove any ~> Centaur metadata expectation metacharacters. */
  private def selectMetadataExpectationSubsetByPrefix(workflow: Workflow, prefix: String, workflowId: WorkflowId, workflowRoot: String): List[(String, JsValue)] = {
    import WorkflowFlatMetadata._
    def replaceVariables(value: JsValue): JsValue = value match {
      case s: JsString => JsString(s.value.replaceExpectationVariables(workflowId, workflowRoot).replaceFirst("^~>", ""))
      case o => o
    }
    val filterLabels: PartialFunction[(String, JsValue), (String, JsValue)] = {
      case (k, v) if k.startsWith(prefix) => k.substring(prefix.length) -> replaceVariables(v)
    }

    for {
      flat <- workflow.metadata.toList
      selected <- flat.value.toList collect filterLabels
    } yield selected
  }

  def validateOutputs(submittedWorkflow: SubmittedWorkflow,
                      workflow: Workflow,
                      workflowRoot: String): Test[Unit] = new Test[Unit] {

    def checkOutputs(expectedOutputs: List[(String, JsValue)])(actualOutputs: Map[String, JsValue]): IO[Unit] = {
      val expected = expectedOutputs.toSet
      val actual = actualOutputs.toSet

      lazy val inActualButNotInExpected = actual.diff(expected)
      lazy val inExpectedButNotInActual = expected.diff(actual)

      if (!workflow.allowOtherOutputs && inActualButNotInExpected.nonEmpty) {
        val message = s"In actual outputs but not in expected and other outputs not allowed: ${inActualButNotInExpected.mkString(", ")}"
        IO.raiseError(CentaurTestException(message, workflow, submittedWorkflow))
      } else if (inExpectedButNotInActual.nonEmpty) {
        val message = s"In expected outputs but not in actual: ${inExpectedButNotInActual.mkString(", ")}"
        IO.raiseError(CentaurTestException(message, workflow, submittedWorkflow))
      } else {
        IO.unit
      }
    }

    override def run: IO[Unit] = {
      import centaur.test.metadata.WorkflowFlatOutputs._

      val expectedOutputs: List[(String, JsValue)] = selectMetadataExpectationSubsetByPrefix(workflow, "outputs.", submittedWorkflow.id, workflowRoot)
      val ioActualOutputs: IO[Map[String, JsValue]] = CentaurCromwellClient.outputs(submittedWorkflow) map { _.asFlat.stringifyValues }

      ioActualOutputs flatMap checkOutputs(expectedOutputs)
    }
  }

  def validateLabels(submittedWorkflow: SubmittedWorkflow,
                     workflow: Workflow,
                     workflowRoot: String): Test[Unit] = new Test[Unit] {
    override def run: IO[Unit] = {
      import centaur.test.metadata.WorkflowFlatLabels._

      val workflowIdLabel = ("cromwell-workflow-id", JsString(s"cromwell-${submittedWorkflow.id}"))
      val expectedLabels: List[(String, JsValue)] = workflowIdLabel ::
        selectMetadataExpectationSubsetByPrefix(workflow, "labels.", submittedWorkflow.id, workflowRoot)
      val ioActualLabels: IO[Map[String, JsValue]] = CentaurCromwellClient.labels(submittedWorkflow) map { _.asFlat.stringifyValues }

      ioActualLabels flatMap { actualLabels =>
        val diff = expectedLabels.toSet.diff(actualLabels.toSet)
        if (diff.nonEmpty) {
          val message = s"In expected labels but not in actual: ${diff.mkString(", ")}"
          IO.raiseError(CentaurTestException(message, workflow, submittedWorkflow))
        } else {
          IO.unit
        }
      }
    }
  }

  /** Compares logs filtered from the raw `metadata` endpoint with the `logs` endpoint. */
  def validateLogs(metadata: WorkflowMetadata, submittedWorkflow: SubmittedWorkflow, workflow: Workflow): Test[Unit] = new Test[Unit] {
    val suffixes = Set("stdout", "shardIndex", "stderr", "attempt", "backendLogs.log")

    def removeSubworkflowKeys(flattened: Map[String, JsValue]): Map[String, JsValue] = {
      val subWorkflowIdPrefixes = flattened.keys.filter(_.endsWith(".subWorkflowId")).map(s => s.substring(0, s.lastIndexOf('.')))
      flattened filter { case (k, _) => !subWorkflowIdPrefixes.exists(k.startsWith) }
    }

    // Filter to only include the fields in the flattened metadata that should appear in the logs endpoint.
    def filterForLogsFields(flattened: Map[String, JsValue]): Map[String, JsValue] = removeSubworkflowKeys(flattened).filter {
      case (k, _) => k == "id" || suffixes.exists(s => k.endsWith("." + s) && !k.contains(".outputs.") && !k.startsWith("outputs."))
    }

    override def run: IO[Unit] = {

      def validateLogsMetadata(flatLogs: Map[String, JsValue], flatFilteredMetadata: Map[String, JsValue]): IO[Unit] =
        if (flatLogs.equals(flatFilteredMetadata)) {
          IO.unit
        } else {
          val message = (List("actual logs endpoint output did not equal filtered metadata", "flat logs: ") ++
            flatLogs.toList ++ List("flat filtered metadata: ") ++ flatFilteredMetadata.toList).mkString("\n")
          IO.raiseError(CentaurTestException(message, workflow, submittedWorkflow))
        }

      for {
        logs <- CentaurCromwellClient.logs(submittedWorkflow)
        flatLogs = logs.asFlat.value
        flatFilteredMetadata = metadata.asFlat.value |> filterForLogsFields
        _ <- validateLogsMetadata(flatLogs, flatFilteredMetadata)
      } yield ()
    }
  }

  def validateMetadata(submittedWorkflow: SubmittedWorkflow,
                       workflowSpec: Workflow,
                       cacheHitUUID: Option[UUID] = None): Test[WorkflowMetadata] = {
    new Test[WorkflowMetadata] {
      def eventuallyMetadata(workflow: SubmittedWorkflow,
                             expectedMetadata: WorkflowFlatMetadata): IO[WorkflowMetadata] = {
        validateMetadata(workflow, expectedMetadata).handleErrorWith({ _ =>
          for {
            _ <- IO.sleep(2.seconds)
            recurse <- eventuallyMetadata(workflow, expectedMetadata)
          } yield recurse
        })
      }

      def validateMetadata(workflow: SubmittedWorkflow,
                           expectedMetadata: WorkflowFlatMetadata): IO[WorkflowMetadata] = {
        def checkDiff(diffs: Iterable[String], actualMetadata: WorkflowMetadata): IO[Unit] = {
          if (diffs.nonEmpty) {
            val message = s"Invalid metadata response:\n -${diffs.mkString("\n -")}\n"
            IO.raiseError(CentaurTestException(message, workflowSpec, workflow, actualMetadata))
          } else {
            IO.unit
          }
        }

        def validateUnwantedMetadata(actualMetadata: WorkflowMetadata): IO[Unit] = {
          if (workflowSpec.notInMetadata.nonEmpty) {
            // Check that none of the "notInMetadata" keys are in the actual metadata
            val absentMdIntersect = workflowSpec.notInMetadata.toSet.intersect(actualMetadata.asFlat.value.keySet)
            if (absentMdIntersect.nonEmpty) {
              val message = s"Found unwanted keys in metadata: ${absentMdIntersect.mkString(", ")}"
              IO.raiseError(CentaurTestException(message, workflowSpec, workflow, actualMetadata))
            } else {
              IO.unit
            }
          } else {
            IO.unit
          }
        }

        def validateAllowOtherOutputs(actualMetadata: WorkflowMetadata): IO[Unit] = {
          if (workflowSpec.allowOtherOutputs) IO.unit
          else {
            val flat = actualMetadata.asFlat.value
            val actualOutputs: Iterable[String] = flat.keys.filter(_.startsWith("outputs."))
            val expectedOutputs: Iterable[String] = workflowSpec.metadata.map(w => w.value.keys.filter(_.startsWith("outputs."))).getOrElse(List.empty)
            val diff = actualOutputs.toSet.diff(expectedOutputs.toSet)
            if (diff.nonEmpty) {
              val message = s"Found unwanted keys in metadata with `allow-other-outputs` = false: ${diff.mkString(", ")}"
              IO.raiseError(CentaurTestException(message, workflowSpec, workflow, actualMetadata))
            } else {
              IO.unit
            }
          }
        }

        for {
          actualMetadata <- CentaurCromwellClient.metadata(workflow)
          _ <- validateUnwantedMetadata(actualMetadata)
          _ <- validateAllowOtherOutputs(actualMetadata)
          diffs = expectedMetadata.diff(actualMetadata.asFlat, workflow.id.id, cacheHitUUID)
          _ <- checkDiff(diffs, actualMetadata)
        } yield actualMetadata
      }

      override def run: IO[WorkflowMetadata] = workflowSpec.metadata match {
        case Some(expectedMetadata) =>
          eventuallyMetadata(submittedWorkflow, expectedMetadata)
            .timeoutTo(CentaurConfig.metadataConsistencyTimeout, validateMetadata(submittedWorkflow, expectedMetadata))
        // Nothing to wait for, so just return the first metadata we get back:
        case None => CentaurCromwellClient.metadata(submittedWorkflow)
      }
    }
  }

  def validateJobManagerStyleMetadata(submittedWorkflow: SubmittedWorkflow,
                                      originalMetadata: String): Test[Unit] = new Test[Unit] {

    def validate(expectation: JsObject, actual: JsObject): IO[Unit] = {
      if(actual.equals(expectation)) {
        IO.pure(())
      } else {

        import diffson._
        import diffson.lcs._
        import diffson.sprayJson._
        import diffson.jsonpatch._
        import diffson.jsonpatch.lcsdiff._

        implicit val lcs = new Patience[JsValue]

        implicit val writer: JsonWriter[JsonPatch[JsValue]] = new JsonWriter[JsonPatch[JsValue]] {
          def processOperation(op: Operation[JsValue]): JsValue = op match {
            case Add(path, value) => JsObject(Map[String, JsValue](
              "op" -> JsString("add"),
              "path" -> JsString(path.toString),
              "value" -> value))
            case Copy(from, path) => JsObject(Map[String, JsValue](
              "op" -> JsString("copy"),
              "from" -> JsString(from.toString),
              "path" -> JsString(path.toString)))
            case Move(from, path) => JsObject(Map[String, JsValue](
              "op" -> JsString("move"),
              "from" -> JsString(from.toString),
              "path" -> JsString(path.toString)))
            case Remove(path, old) => JsObject(Map[String, JsValue](
              "op" -> JsString("remove"),
              "path" -> JsString(path.toString)) ++ old.map(o => "old" -> o) )
            case Replace(path, value, old) => JsObject(Map[String, JsValue](
              "op" -> JsString("replace"),
              "path" -> JsString(path.toString),
              "value" -> value) ++ old.map(o => "old" -> o) )
            case diffson.jsonpatch.Test(path, value) => JsObject(Map[String, JsValue](
              "op" -> JsString("test"),
              "path" -> JsString(path.toString),
              "value" -> value))
          }

          override def write(obj: JsonPatch[JsValue]): JsValue = {
            JsArray(obj.ops.toVector.map(processOperation))
          }
        }

        val jsonDiff = diff(expectation: JsValue, actual: JsValue)

        logger.error(s"Bad JM style metadata: ${jsonDiff.toJson.prettyPrint}")
        IO.raiseError(new Exception(s"Bad JM style metadata. See error log output"))
      }
    }

    override def run: IO[Unit] = for {
      jmMetadata <- CentaurCromwellClient.metadata(workflow = submittedWorkflow, Option(CentaurCromwellClient.defaultMetadataArgs.getOrElse(Map.empty) ++ jmArgs))
      jmMetadataObject <- IO.fromTry(Try(jmMetadata.value.parseJson.asJsObject))
      expectation <- IO.fromTry(Try(extractJmStyleMetadataFields(originalMetadata.parseJson.asJsObject)))

      validationUnit <- validate(expectation, jmMetadataObject)
    } yield validationUnit
  }

  val oneWordIncludeKeys = List(
    "attempt", "callRoot", "end",
    "executionStatus", "failures", "inputs", "jobId",
    "calls", "outputs", "shardIndex", "start", "stderr", "stdout",
    "description", "executionEvents", "labels", "parentWorkflowId",
    "returnCode", "status", "submission", "subWorkflowId", "workflowName"
  )

  val jmArgs = Map(
    "includeKey" -> (oneWordIncludeKeys :+ "callCaching:hit"),
    "excludeKey" -> List("callCaching:hitFailures"),
    "expandSubWorkflows" -> List("false")
  )

  // I chose to re-implement "process metadata for job manager's includeKeys/excludeKeys" and not re-use production.
  // Reasoning is:
  //  1. the ultra-specific JM case is a lot simpler to reason about and can be done in a few lines
  //  2. it's nice to be able to assert general cases rather than singular specific examples
  //  3. it's hopefully unlikely that the same bug will show up in two separate implementations
  def extractJmStyleMetadataFields(originalWorkflowMetadataJson: JsObject): JsObject = {
    val originalCallMetadataJson = originalWorkflowMetadataJson.fields.get("calls").map(_.asJsObject)

    // NB: this filter to remove "calls" is because - although it is a single word in the JM request,
    // it gets treated specially by the API (so has to be treated specially here too)
    def processOneWordIncludes(json: JsObject) = (oneWordIncludeKeys.filterNot(_ == "calls") :+ "id").foldRight(JsObject.empty) { (toInclude, current) =>
      json.fields.get(toInclude) match {
        case Some(jsonToInclude) => JsObject(current.fields + (toInclude -> jsonToInclude))
        case None => current
      }
    }

    def processCallCacheField(callJson: JsObject) = for {
      originalCallCachingField <- callJson.fields.get("callCaching")
      originalHitField <- originalCallCachingField.asJsObject.fields.get("hit")
    } yield "callCaching" -> JsObject(Map("hit" -> originalHitField))

    def processCallsSection(calls: JsObject): JsObject = {
      val newFields = calls.fields.map { case (callName, callsArray) =>
        val newElements = callsArray.asInstanceOf[spray.json.JsArray].elements.map(_.asJsObject).map { call =>
          val callWithOneWordIncludes = processOneWordIncludes(call)
          val callCacheField = processCallCacheField(call)
          JsObject(callWithOneWordIncludes.fields ++ callCacheField)
        }
        callName -> JsArray(newElements)
      }
      JsObject(newFields)
    }

    val workflowLevelWithOneWordIncludes = processOneWordIncludes(originalWorkflowMetadataJson)
    val callsField = originalCallMetadataJson map { calls => Map("calls" -> processCallsSection(calls)) } getOrElse Map.empty

    JsObject(workflowLevelWithOneWordIncludes.fields ++ callsField)
  }

  def waitForArchive(submittedWorkflow: SubmittedWorkflow, workflowDefinition: Workflow): Test[Unit] = {
    new Test[Unit] {

      def validateMetadataArchiveStatus(status: String): IO[Boolean] = {
        logger.info(s"Validating archive status '$status for workflow ID: ${submittedWorkflow.id}'")
        if (status == "Archived") {
          IO.pure(true)
        }  else if (status == "Unarchived" ) {
          IO.pure(false)
        } else {
          IO.fromTry(Failure(new Exception(s"Expected Archived but got $status")))
        }
      }

      def checkArchived(): IO[(Boolean, Boolean)] = for {
        archiveStatus <- CentaurCromwellClient.archiveStatus(submittedWorkflow.id)
        isArchived <- validateMetadataArchiveStatus(archiveStatus)
        isMetadataSourceArchived <- validateMetadataSourceArchived()
      } yield (isArchived, isMetadataSourceArchived)

      def validateMetadataSourceArchived(): IO[Boolean] = for {
        metadataSource <- CentaurCromwellClient.metadataWithId(submittedWorkflow.id).map(_.asFlat.stringifyValues.get("metadataSource"))
        isMetadataSourceArchived <- {
          if (metadataSource.contains(JsString("Archived"))) {
            IO.pure(true)
          } else if (metadataSource.contains(JsString("Unarchived"))) {
            IO.pure(false)
          } else {
            IO.raiseError(CentaurTestException(
              s"`Metadata` endpoint returned unknown value for `metadataSource`: $metadataSource",
              workflowDefinition,
              submittedWorkflow
            ))
          }
        }
      } yield isMetadataSourceArchived

      def eventuallyArchived(): IO[Unit] = {
        checkArchived() flatMap {
          case (true, true) => IO.pure(())
          case (true, false) =>
            IO.raiseError(CentaurTestException(
              "`Query` endpoint returns metadata status \"Archived\" but `metadata` endpoint returns metadata source \"Unarchived\"",
              workflowDefinition,
              submittedWorkflow
            ))
          case (false, true) | (false, false) => for {
            _ <- IO.sleep(2.seconds)
            recurse <- eventuallyArchived()
          } yield recurse
        }
      }

      override def run: IO[Unit] = {
        if (CentaurConfig.expectCarbonite) {
          eventuallyArchived().timeout(CentaurConfig.metadataConsistencyTimeout)
        } else {
          IO.pure(())
        }
      }
    }
  }

  /**
    * Verify that none of the calls within the workflow are cached.
    */
  def validateCacheResultField(workflowDefinition: Workflow,
                               submittedWorkflow: SubmittedWorkflow,
                               metadata: WorkflowMetadata,
                               blacklistedValue: String): Test[Unit] = {
    new Test[Unit] {
      override def run: IO[Unit] = {
        val badCacheResults = metadata.asFlat.value collect {
          case (k, JsString(v)) if k.contains("callCaching.result") && v.contains(blacklistedValue) => s"$k: $v"
        }

        if (badCacheResults.isEmpty) {
          IO.unit
        } else {
          val message = s"Found unexpected cache hits for " +
            s"${workflowDefinition.testName}:${badCacheResults.mkString("\n", "\n", "\n")}"
          IO.raiseError(CentaurTestException(message, workflowDefinition, submittedWorkflow, metadata))
        }
      }
    }
  }

  def validateDirectoryContentsCounts(workflowDefinition: Workflow,
                                      submittedWorkflow: SubmittedWorkflow,
                                      metadata: WorkflowMetadata): Test[Unit] = new Test[Unit] {
    private val workflowId = submittedWorkflow.id.id.toString

    override def run: IO[Unit] = workflowDefinition.directoryContentCounts match {
      case None => IO.unit
      case Some(directoryContentCountCheck) =>
        val counts = directoryContentCountCheck.expectedDirectoryContentsCounts map {
          case (directory, count) =>
            val substitutedDir = directory.replaceAll("<<UUID>>", workflowId)
            (substitutedDir, count, directoryContentCountCheck.checkFiles.countObjectsAtPath(substitutedDir))
        }

        val badCounts = counts collect {
          case (directory, expectedCount, actualCount) if expectedCount != actualCount => s"Expected to find $expectedCount item(s) at $directory but got $actualCount"
        }
        if (badCounts.isEmpty) {
          IO.unit
        } else {
          val message = badCounts.mkString("\n", "\n", "\n")
          IO.raiseError(CentaurTestException(message, workflowDefinition, submittedWorkflow, metadata))
        }
    }
  }

  def validateNoCacheHits(submittedWorkflow: SubmittedWorkflow,
                          metadata: WorkflowMetadata,
                          workflowDefinition: Workflow): Test[Unit] = {
    validateCacheResultField(workflowDefinition, submittedWorkflow, metadata, "Cache Hit")
  }

  def validateNoCacheMisses(submittedWorkflow: SubmittedWorkflow,
                            metadata: WorkflowMetadata,
                            workflowDefinition: Workflow): Test[Unit] = {
    validateCacheResultField(workflowDefinition, submittedWorkflow, metadata, "Cache Miss")
  }

  def validateSubmitFailure(workflow: Workflow,
                            expectedSubmitResponse: SubmitHttpResponse,
                            actualSubmitResponse: SubmitHttpResponse): Test[Unit] = {
    new Test[Unit] {
      override def run: IO[Unit] = {
        if (expectedSubmitResponse == actualSubmitResponse) {
          IO.unit
        } else {
          val message =
            s"""|
                |Expected
                |$expectedSubmitResponse
                |
                |but got:
                |$actualSubmitResponse
                |""".stripMargin
          IO.raiseError(CentaurTestException(message, workflow))
        }
      }
    }
  }
}
