package centaur

import better.files._
import cats.data.Validated.{Invalid, Valid}
import cats.effect.IO
import cats.instances.list._
import cats.syntax.flatMap._
import cats.syntax.functor._
import cats.syntax.traverse._
import centaur.reporting.{ErrorReporters, TestEnvironment}
import centaur.test.CentaurTestException
import centaur.test.standard.CentaurTestCase
import org.scalatest._

import scala.concurrent.Future

@DoNotDiscover
abstract class AbstractCentaurTestCaseSpec(cromwellBackends: List[String]) extends AsyncFlatSpec with Matchers {

  /*
  NOTE: We need to statically initialize the object so that the exceptions appear here in the class constructor.
  Otherwise instead of seeing 'java.lang.ExceptionInInitializerError' with a stack trace, there will be no tests
  generated, and no errors, resulting in 'Total number of tests run: 0' and 'No tests were executed.'
   */
  ErrorReporters.getClass

  private def testCases(baseFile: File): List[CentaurTestCase] = {
    val files = baseFile.list.filter(_.isRegularFile).toList
    val testCases = files.traverse(CentaurTestCase.fromFile)

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

    def runTest(): IO[Unit] = testCase.testFunction.run.void

    // Make tags, but enforce lowercase:
    val tags = (testCase.testOptions.tags :+ testCase.workflow.testName :+ testCase.testFormat.name) map { x => Tag(x.toLowerCase) }
    val isIgnored = testCase.isIgnored(cromwellBackends)

    runOrDont(nameTest, tags, isIgnored, runTest())
  }

  def executeUpgradeTest(testCase: CentaurTestCase): Unit =
    executeStandardTest(upgradeTestWdl(testCase))

  private def upgradeTestWdl(testCase: CentaurTestCase): CentaurTestCase = {
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
          zippedImports = Option(upgradedImportsDir.zip())))) // An empty zip appears to be completely harmless, so no special handling

    rootWorkflowFile.delete(true)
    upgradedImportsDir.delete(true)
    workingDir.delete(true)

    newCase
  }

  private def runOrDont(testName: String, tags: List[Tag], ignore: Boolean, runTest: => IO[Unit]): Unit = {

    val itShould: ItVerbString = it should testName

    tags match {
      case Nil => runOrDont(itShould, ignore, testName, runTest)
      case head :: Nil => runOrDont(itShould taggedAs head, ignore, testName, runTest)
      case head :: tail => runOrDont(itShould taggedAs(head, tail: _*), ignore, testName, runTest)
    }
  }

  private def runOrDont(itVerbString: ItVerbString,
                        ignore: Boolean,
                        testName: String,
                        runTest: => IO[Unit]): Unit = {
    if (ignore) {
      itVerbString ignore Future.successful(succeed)
    } else {
      itVerbString in tryTryAgain(testName, runTest, ErrorReporters.retryAttempts).unsafeToFuture().map(_ => succeed)
    }
  }

  private def runOrDont(itVerbStringTaggedAs: ItVerbStringTaggedAs,
                        ignore: Boolean,
                        testName: String,
                        runTest: => IO[Unit]): Unit = {
    if (ignore) {
      itVerbStringTaggedAs ignore Future.successful(succeed)
    } else {
      itVerbStringTaggedAs in
        tryTryAgain(testName, runTest, ErrorReporters.retryAttempts).unsafeToFuture().map(_ => succeed)
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
  private def tryTryAgain(testName: String, runTest: => IO[Unit], retries: Int, attempt: Int = 0): IO[Unit] = {
    def maybeRetry(centaurTestException: CentaurTestException): IO[Unit] = {
      val testEnvironment = TestEnvironment(testName, retries, attempt)
      for {
        _ <- ErrorReporters.logCentaurFailure(testEnvironment, centaurTestException)
        _ <- if (attempt < retries) {
          tryTryAgain(testName, runTest, retries, attempt + 1)
        } else {
          IO.raiseError(centaurTestException)
        }
      } yield ()
    }

    val runTestIo = IO(runTest).flatten

    runTestIo handleErrorWith {
      case centaurTestException: CentaurTestException => maybeRetry(centaurTestException)
      case _ => runTestIo
    }
  }

}
