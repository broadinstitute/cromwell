package cromwell

import akka.testkit._
import cromwell.util.SampleWdl
import wom.values.WomString


class PostfixQuantifierWorkflowSpec extends CromwellTestKitWordSpec {
  "A task which contains a parameter with a zero-or-more postfix quantifier" should {
    "accept an array of size 3" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ZeroOrMorePostfixQuantifierWorkflowWithArrayInput,
        EventFilter.info(pattern = "Starting calls: postfix.hello", occurrences = 1),
        expectedOutputs = Map("postfix.hello.greeting" -> WomString("hello alice,bob,charles"))
      )
    }
    "accept an array of size 1" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ZeroOrMorePostfixQuantifierWorkflowWithOneElementArrayInput,
        EventFilter.info(pattern = "Starting calls: postfix.hello", occurrences = 1),
        expectedOutputs = Map("postfix.hello.greeting" -> WomString("hello alice"))
      )
    }
    "accept an array of size 0" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ZeroOrMorePostfixQuantifierWorkflowWithZeroElementArrayInput,
        EventFilter.info(pattern = "Starting calls: postfix.hello", occurrences = 1),
        expectedOutputs = Map("postfix.hello.greeting" -> WomString("hello"))
      )
    }
  }

  "A task which contains a parameter with a one-or-more postfix quantifier" should {
    "accept an array for the value" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.OneOrMorePostfixQuantifierWorkflowWithArrayInput,
        EventFilter.info(pattern = "Starting calls: postfix.hello", occurrences = 1),
        expectedOutputs = Map("postfix.hello.greeting" -> WomString("hello alice,bob,charles"))
      )
    }
    "accept a scalar for the value" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.OneOrMorePostfixQuantifierWorkflowWithScalarInput,
        EventFilter.info(pattern = "Starting calls: postfix.hello", occurrences = 1),
        expectedOutputs = Map("postfix.hello.greeting" -> WomString("hello alice"))
      )
    }
  }
}
