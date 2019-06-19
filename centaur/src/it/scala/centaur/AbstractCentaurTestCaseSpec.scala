package centaur

import better.files._
import cats.data.Validated.{Invalid, Valid}
import cats.effect.IO
import cats.instances.list._
import cats.syntax.flatMap._
import cats.syntax.traverse._
import centaur.reporting.{ErrorReporters, SuccessReporters, TestEnvironment}
import centaur.test.CentaurTestException
import centaur.test.standard.CentaurTestCase
import centaur.test.submit.{SubmitResponse, SubmitWorkflowResponse}
import org.scalatest._

import scala.concurrent.Future

@DoNotDiscover
abstract class AbstractCentaurTestCaseSpec(cromwellBackends: List[String], cromwellTracker: Option[CromwellTracker] = None) extends AsyncFlatSpec with Matchers {

  /*
  NOTE: We need to statically initialize the object so that the exceptions appear here in the class constructor.
  Otherwise instead of seeing 'java.lang.ExceptionInInitializerError' with a stack trace, there will be no tests
  generated, and no errors, resulting in 'Total number of tests run: 0' and 'No tests were executed.'
   */
  ErrorReporters.getClass
  SuccessReporters.getClass

  private def testCases(baseFile: File): List[CentaurTestCase] = {
    val files = baseFile.list.filter(_.isRegularFile).toList
    val testCases = files.traverse(CentaurTestCase.fromFile(cromwellTracker))

    testCases match {
      case Valid(l) => l
      case Invalid(e) => throw new IllegalStateException("\n" + e.toList.mkString("\n") + "\n")
    }
  }

  def allTestCases: List[CentaurTestCase] = {
    val optionalTestCases = CentaurConfig.optionalTestPath map (File(_)) map testCases getOrElse List.empty
    val standardTestCases = testCases(CentaurConfig.standardTestCasePath)
    optionalTestCases ++ standardTestCases
  }

  def executeStandardTest(testCase: CentaurTestCase): Unit = {
    def nameTest = s"${testCase.testFormat.testSpecString} ${testCase.workflow.testName}"

    def runTest(): IO[SubmitResponse] = testCase.testFunction.run

    // Make tags, but enforce lowercase:
    val tags = (testCase.testOptions.tags :+ testCase.workflow.testName :+ testCase.testFormat.name) map { x => Tag(x.toLowerCase) }
    val isIgnored = testCase.isIgnored(cromwellBackends)
    val retries = if (testCase.workflow.retryTestFailures) ErrorReporters.retryAttempts else 0

    runOrDont(nameTest, tags, isIgnored, retries, runTest())
  }

  def executeWdlUpgradeTest(testCase: CentaurTestCase): Unit =
    executeStandardTest(wdlUpgradeTestWdl(testCase))

  private def wdlUpgradeTestWdl(testCase: CentaurTestCase): CentaurTestCase = {
    import better.files.File
    import womtool.WomtoolMain

    // The suffix matters because WomGraphMaker.getBundle() uses it to choose the language factory
    val rootWorkflowFile = File.newTemporaryFile(suffix = ".wdl").append(testCase.workflow.data.workflowContent.get)
    val workingDir = File.newTemporaryDirectory()
    val upgradedImportsDir = File.newTemporaryDirectory()
    val rootWorkflowFilepath = workingDir / rootWorkflowFile.name

    // Un-upgraded imports go into the working directory
    testCase.workflow.data.zippedImports match {
      case Some(importsZip: File) =>
        importsZip.unzipTo(workingDir)
      case None => ()
    }

    // Upgrade the imports and copy to main working dir (precludes transitive imports; no recursion yet)
    workingDir.list.toList.map { file: File =>
      val upgradedWdl = WomtoolMain.upgrade(file.pathAsString).stdout.get
      upgradedImportsDir.createChild(file.name).append(upgradedWdl)
    }

    // Copy to working directory after we operate on the imports that are in it
    rootWorkflowFile.copyTo(rootWorkflowFilepath)

    val upgradeResult = WomtoolMain.upgrade(rootWorkflowFilepath.pathAsString)

    upgradeResult.stderr match {
      case Some(stderr) => println(stderr)
      case _ => ()
    }

    val newCase = testCase.copy(
      workflow = testCase.workflow.copy(
        testName = testCase.workflow.testName + " (draft-2 to 1.0 upgrade)",
        data = testCase.workflow.data.copy(
          workflowContent = Option(upgradeResult.stdout.get), // this '.get' catches an error if upgrade fails
          zippedImports = Option(upgradedImportsDir.zip()))))(cromwellTracker) // An empty zip appears to be completely harmless, so no special handling

    rootWorkflowFile.delete(true)
    upgradedImportsDir.delete(true)
    workingDir.delete(true)

    newCase
  }

  private def runOrDont(testName: String,
                        tags: List[Tag],
                        ignore: Boolean,
                        retries: Int,
                        runTest: => IO[SubmitResponse]): Unit = {

    val itShould: ItVerbString = it should testName

    tags match {
      case Nil => runOrDont(itShould, ignore, testName, retries, runTest)
      case head :: Nil => runOrDont(itShould taggedAs head, ignore, testName, retries, runTest)
      case head :: tail => runOrDont(itShould taggedAs(head, tail: _*), ignore, testName, retries, runTest)
    }
  }

  private def runOrDont(itVerbString: ItVerbString,
                        ignore: Boolean,
                        testName: String,
                        retries: Int,
                        runTest: => IO[SubmitResponse]): Unit = {
    if (ignore) {
      itVerbString ignore Future.successful(succeed)
    } else {
      itVerbString in tryTryAgain(testName, runTest, retries).unsafeToFuture().map(_ => succeed)
    }
  }

  private def runOrDont(itVerbStringTaggedAs: ItVerbStringTaggedAs,
                        ignore: Boolean,
                        testName: String,
                        retries: Int,
                        runTest: => IO[SubmitResponse]): Unit = {
    if (ignore) {
      itVerbStringTaggedAs ignore Future.successful(succeed)
    } else {
      itVerbStringTaggedAs in
        tryTryAgain(testName, runTest, retries).unsafeToFuture().map(_ => succeed)
    }
  }

  /**
    * Returns an IO effect that will recursively try to run a test.
    *
    * @param testName Name of the ScalaTest.
    * @param runTest Thunk to run the test.
    * @param retries Total number of attempts to retry.
    * @param attempt Current zero based attempt.
    * @return IO effect that will run the test, possibly retrying.
    */
  private def tryTryAgain(testName: String, runTest: => IO[SubmitResponse], retries: Int, attempt: Int = 0): IO[SubmitResponse] = {
    def maybeRetry(centaurTestException: CentaurTestException): IO[SubmitResponse] = {
      val testEnvironment = TestEnvironment(testName, retries, attempt)
      for {
        _ <- ErrorReporters.logCentaurFailure(testEnvironment, centaurTestException)
        r <- if (attempt < retries) {
          tryTryAgain(testName, runTest, retries, attempt + 1)
        } else {
          IO.raiseError(centaurTestException)
        }
      } yield r
    }

    val runTestIo = IO(runTest).flatten

    runTestIo.redeemWith(
      {
        case centaurTestException: CentaurTestException => maybeRetry(centaurTestException)
        case _ => runTestIo
      },
      {
        case workflowResponse: SubmitWorkflowResponse => SuccessReporters.logSuccessfulRun(workflowResponse)
        case other => IO.pure(other)
      }
    )
  }

}
