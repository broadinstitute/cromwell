package centaur

import java.nio.file.Path

import cats.data.Validated.{Invalid, Valid}
import cats.instances.list._
import cats.syntax.traverse._
import centaur.test.standard.CentaurTestCase
import org.scalatest.{DoNotDiscover, FlatSpec, Matchers, Tag}
import wdl.draft2.model.WdlNamespace
import wdl.transforms.draft2.wdlom2wom.WdlDraft2WomBundleMakers.wdlDraft2NamespaceWomBundleMaker

@DoNotDiscover
abstract class AbstractCentaurTestCaseSpec(cromwellBackends: List[String]) extends FlatSpec with Matchers {

  private def testCases(basePath: Path): List[CentaurTestCase] = {
    val files = basePath.toFile.listFiles.toList collect { case x if x.isFile => x.toPath }
    val testCases = files.traverse(CentaurTestCase.fromPath)

    testCases match {
      case Valid(l) => l
      case Invalid(e) => throw new IllegalStateException("\n" + e.toList.mkString("\n") + "\n")
    }
  }
  
  def allTestCases: List[CentaurTestCase] = {
    val optionalTestCases = CentaurConfig.optionalTestPath map testCases getOrElse List.empty
    val standardTestCases = testCases(CentaurConfig.standardTestCasePath)
    optionalTestCases ++ standardTestCases
  }

  def executeStandardTest(testCase: CentaurTestCase): Unit = {
    def nameTest = s"${testCase.testFormat.testSpecString} ${testCase.workflow.testName}"
//    def runTest(): Unit = {
//      testCase.testFunction.run.get
//      ()
//    }

    // Make tags, but enforce lowercase:
    val tags = (testCase.testOptions.tags :+ testCase.workflow.testName :+ testCase.testFormat.name) map { x => Tag(x.toLowerCase) }
    val isIgnored = testCase.isIgnored(cromwellBackends)

//    if (!(testCase.workflow.data.workflowContent.contains("version 1.0") ||
//          testCase.workflow.data.workflowContent.contains("version draft-3")) &&
//        !testCase.workflow.data.workflowType.contains("CWL") &&
//        !tags.map(_.name).contains("subworkflow")) {
    if (tags.map(_.name).contains("upgrade")) {
      val draft2Wdl = WdlNamespace.loadUsingSource(
        testCase.workflow.data.workflowContent,
        None,
        Some(Seq(WdlNamespace.fileResolver))
      ).get

      val bundle = wdlDraft2NamespaceWomBundleMaker.toWomBundle(draft2Wdl).right.get

      import wdl.draft3.transforms.wdlom2wdl.WdlWriter.ops._
      import wdl.draft3.transforms.wdlom2wdl.WdlWriterImpl.fileElementWriter
      import wdl.draft3.transforms.wom2wdlom.WomToWdlom.ops._
      import wdl.draft3.transforms.wom2wdlom.WomToWdlomImpl.womBundleToFileElement

      val v1testCase = testCase.copy(
        workflow = testCase.workflow.copy(
          data = testCase.workflow.data.copy(workflowContent = bundle.toWdlom.toWdlV1)))

      def runTest(): Unit = {
        v1testCase.testFunction.run.get
        ()
      }

      runOrDont(nameTest + "_upgrade", tags, isIgnored, runTest())
    }

//    runOrDont(nameTest, tags, isIgnored, runTest())
  }

  private def runOrDont(testName: String, tags: List[Tag], ignore: Boolean, runTest: => Any): Unit = {

    val itShould: ItVerbString = it should testName

    tags match {
      case Nil => runOrDont(itShould, ignore, runTest)
      case head :: Nil => runOrDont(itShould taggedAs head, ignore, runTest)
      case head :: tail => runOrDont(itShould taggedAs(head, tail: _*), ignore, runTest)
    }
  }

  private def runOrDont(itVerbString: ItVerbString, ignore: Boolean, runTest: => Any): Unit = {
    if (ignore) {
      itVerbString ignore runTest
    } else {
      itVerbString in runTest
    }
  }

  private def runOrDont(itVerbStringTaggedAs: ItVerbStringTaggedAs, ignore: Boolean, runTest: => Any): Unit = {
    if (ignore) {
      itVerbStringTaggedAs ignore runTest
    } else {
      itVerbStringTaggedAs in runTest
    }
  }
  
}
