package cromwell

import akka.testkit._
import wdl4s.values.WdlString
import cromwell.util.SampleWdl


class PostfixQuantifierWorkflowSpec extends CromwellTestKitSpec {
  "A task which contains a parameter with a zero-or-more postfix quantifier" should {
    "accept an array of size 3" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ZeroOrMorePostfixQuantifierWorkflowWithArrayInput,
        EventFilter.info(pattern = "Starting calls: postfix.hello", occurrences = 1),
        expectedOutputs = Map("postfix_hello_greeting" -> WdlString("hello alice,bob,charles"))
      )
    }
    "accept an array of size 1" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ZeroOrMorePostfixQuantifierWorkflowWithOneElementArrayInput,
        EventFilter.info(pattern = "Starting calls: postfix.hello", occurrences = 1),
        expectedOutputs = Map("postfix_hello_greeting" -> WdlString("hello alice"))
      )
    }
    "accept an array of size 0" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ZeroOrMorePostfixQuantifierWorkflowWithZeroElementArrayInput,
        EventFilter.info(pattern = "Starting calls: postfix.hello", occurrences = 1),
        expectedOutputs = Map("postfix_hello_greeting" -> WdlString("hello"))
      )
    }
  }

  "A task which contains a parameter with a one-or-more postfix quantifier" should {
    "accept an array for the value" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.OneOrMorePostfixQuantifierWorkflowWithArrayInput,
        EventFilter.info(pattern = "Starting calls: postfix.hello", occurrences = 1),
        expectedOutputs = Map("postfix_hello_greeting" -> WdlString("hello alice,bob,charles"))
      )
    }
    "accept a scalar for the value" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.OneOrMorePostfixQuantifierWorkflowWithScalarInput,
        EventFilter.info(pattern = "Starting calls: postfix.hello", occurrences = 1),
        expectedOutputs = Map("postfix_hello_greeting" -> WdlString("hello alice"))
      )
    }
  }
}
