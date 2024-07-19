package centaur

import better.files._
import cats.data.Validated.{Invalid, Valid}
import cats.effect.IO
import cats.instances.list._
import cats.syntax.flatMap._
import cats.syntax.traverse._
import centaur.callcaching.CromwellDatabaseCallCaching
import centaur.reporting.{ErrorReporters, SuccessReporters, TestEnvironment}
import centaur.test.CentaurTestException
import centaur.test.standard.CentaurTestCase
import centaur.test.submit.{SubmitResponse, SubmitWorkflowResponse}
import centaur.test.workflow.WorkflowData
import cromwell.api.model.WorkflowId
import org.scalatest._
import org.scalatest.flatspec.AsyncFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.concurrent.Future

@DoNotDiscover
abstract class AbstractCentaurTestCaseSpec(cromwellBackends: List[String],
                                           cromwellTracker: Option[CromwellTracker] = None
) extends AsyncFlatSpec
    with Matchers {

  /*
  NOTE: We need to statically initialize the object so that the exceptions appear here in the class constructor.
  Otherwise instead of seeing 'java.lang.ExceptionInInitializerError' with a stack trace, there will be no tests
  generated, and no errors, resulting in 'Total number of tests run: 0' and 'No tests were executed.'
   */
  ErrorReporters.getClass
  SuccessReporters.getClass

  private def testCases(baseFile: File): List[CentaurTestCase] = {
    val files = baseFile.listRecursively.filter(isTestFile).toList
    val testCases = files.traverse(CentaurTestCase.fromFile(cromwellTracker))

    testCases match {
      case Valid(l) => l
      case Invalid(e) => throw new IllegalStateException("\n" + e.toList.mkString("\n") + "\n")
    }
  }

  private def isTestFile(file: File) =
    file.isRegularFile && file.extension.contains(".test")

  def allTestCases: List[CentaurTestCase] = {
    val optionalTestCases = CentaurConfig.optionalTestPath map (File(_)) map testCases getOrElse List.empty
    val standardTestCases = testCases(CentaurConfig.standardTestCasePath)
    val allTestsCases = optionalTestCases ++ standardTestCases
    val duplicateTestNames = allTestsCases
      .map(_.workflow.testName)
      .groupBy(identity)
      .collect { case (key, values) if values.lengthCompare(1) > 0 => key }
    if (duplicateTestNames.nonEmpty) {
      throw new RuntimeException(
        "The following test names are duplicated in more than one test file: " +
          duplicateTestNames.mkString(", ")
      )
    }
    allTestsCases
  }

  def executeStandardTest(testCase: CentaurTestCase): Unit = {
    def runTestAndDeleteZippedImports(): IO[SubmitResponse] = for {
      submitResponse <- testCase.testFunction.run
      _ = cleanUpImports(testCase.workflow.data) // cleanup imports in the end of test
    } yield submitResponse

    // Make tags, but enforce lowercase:
    val tags = (testCase.testOptions.tags :+ testCase.workflow.testName :+ testCase.testFormat.name) map { x =>
      Tag(x.toLowerCase)
    }
    val isIgnored = testCase.isIgnored(cromwellBackends)
    val retries = if (testCase.workflow.retryTestFailures) ErrorReporters.retryAttempts else 0

    runOrDont(testCase, tags, isIgnored, retries, runTestAndDeleteZippedImports())
  }

  private def runOrDont(testCase: CentaurTestCase,
                        tags: List[Tag],
                        ignore: Boolean,
                        retries: Int,
                        runTest: => IO[SubmitResponse]
  ): Unit = {

    val itShould: ItVerbString = it should testCase.name

    tags match {
      case Nil => runOrDont(itShould, ignore, testCase, retries, runTest)
      case head :: Nil => runOrDont(itShould taggedAs head, ignore, testCase, retries, runTest)
      case head :: tail => runOrDont(itShould taggedAs (head, tail: _*), ignore, testCase, retries, runTest)
    }
  }

  private def runOrDont(itVerbString: ItVerbString,
                        ignore: Boolean,
                        testCase: CentaurTestCase,
                        retries: Int,
                        runTest: => IO[SubmitResponse]
  ): Unit =
    if (ignore) {
      itVerbString ignore Future.successful(succeed)
    } else {
      itVerbString in tryTryAgain(testCase, runTest, retries).unsafeToFuture().map(_ => succeed)
    }

  private def runOrDont(itVerbStringTaggedAs: ItVerbStringTaggedAs,
                        ignore: Boolean,
                        testCase: CentaurTestCase,
                        retries: Int,
                        runTest: => IO[SubmitResponse]
  ): Unit =
    if (ignore) {
      itVerbStringTaggedAs ignore Future.successful(succeed)
    } else {
      itVerbStringTaggedAs in
        tryTryAgain(testCase, runTest, retries).unsafeToFuture().map(_ => succeed)
    }

  /**
    * Returns an IO effect that will recursively try to run a test.
    *
    * @param testCase CentaurTestCase.
    * @param runTest Thunk to run the test.
    * @param retries Total number of attempts to retry.
    * @param attempt Current zero based attempt.
    * @return IO effect that will run the test, possibly retrying.
    */
  private def tryTryAgain(testCase: CentaurTestCase,
                          runTest: => IO[SubmitResponse],
                          retries: Int,
                          attempt: Int = 0
  ): IO[SubmitResponse] = {

    def maybeRetry(centaurTestException: CentaurTestException): IO[SubmitResponse] = {

      def clearCachedResults(workflowId: WorkflowId): IO[Unit] =
        CromwellDatabaseCallCaching.clearCachedResults(workflowId.toString)

      val testEnvironment = TestEnvironment(testCase, retries, attempt)
      for {
        _ <- ErrorReporters.logFailure(testEnvironment, centaurTestException)
        r <-
          if (attempt < retries) {
            testCase.submittedWorkflowTracker.cleanUpBeforeRetry(clearCachedResults) *>
              tryTryAgain(testCase, runTest, retries, attempt + 1)
          } else {
            IO.raiseError(centaurTestException)
          }
      } yield r
    }

    val runTestIo = IO(runTest).flatten

    runTestIo.redeemWith(
      {
        case centaurTestException: CentaurTestException =>
          maybeRetry(centaurTestException)
        case nonCentaurThrowable: Throwable =>
          val testEnvironment = TestEnvironment(testCase, retries = attempt + 1, attempt) // allow one last retry
          ErrorReporters.logFailure(testEnvironment, nonCentaurThrowable)
          runTestIo
      },
      {
        case workflowResponse: SubmitWorkflowResponse => SuccessReporters.logSuccessfulRun(workflowResponse)
        case other => IO.pure(other)
      }
    )
  }

  /**
   * Clean up temporary zip files created for Imports testing.
   */
  private def cleanUpImports(wfData: WorkflowData) =
    wfData.zippedImports match {
      case Some(zipFile) => zipFile.delete(swallowIOExceptions = true)
      case None => //
    }

}
