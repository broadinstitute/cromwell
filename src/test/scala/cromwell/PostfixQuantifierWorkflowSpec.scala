package cromwell

import akka.testkit._
import cromwell.binding.values.WdlString
import cromwell.util.SampleWdl

import scala.language.postfixOps

class PostfixQuantifierWorkflowSpec extends CromwellTestkitSpec("PostfixQuantifierWorkflowSpec") {
  "A task which contains a parameter with a zero-or-more postfix quantifier" should {
    "accept an array of size 3" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ZeroOrMorePostfixQuantifierWorkflowWithArrayInput,
        EventFilter.info(pattern = s"starting calls: postfix.hello", occurrences = 1),
        expectedOutputs = Map("postfix.hello.greeting" -> WdlString("hello alice,bob,charles"))
      )
    }
    "accept an array of size 1" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ZeroOrMorePostfixQuantifierWorkflowWithOneElementArrayInput,
        EventFilter.info(pattern = s"starting calls: postfix.hello", occurrences = 1),
        expectedOutputs = Map("postfix.hello.greeting" -> WdlString("hello alice"))
      )
    }
    "accept an array of size 0" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ZeroOrMorePostfixQuantifierWorkflowWithZeroElementArrayInput,
        EventFilter.info(pattern = s"starting calls: postfix.hello", occurrences = 1),
        expectedOutputs = Map("postfix.hello.greeting" -> WdlString("hello "))
      )
    }
  }

  "A task which contains a parameter with a one-or-more postfix quantifier" should {
    "accept an array for the value" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.OneOrMorePostfixQuantifierWorkflowWithArrayInput,
        EventFilter.info(pattern = s"starting calls: postfix.hello", occurrences = 1),
        expectedOutputs = Map("postfix.hello.greeting" -> WdlString("hello alice,bob,charles"))
      )
    }
    "accept a scalar for the value" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.OneOrMorePostfixQuantifierWorkflowWithScalarInput,
        EventFilter.info(pattern = s"starting calls: postfix.hello", occurrences = 1),
        expectedOutputs = Map("postfix.hello.greeting" -> WdlString("hello alice"))
      )
    }
  }
}
