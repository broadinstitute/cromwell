package wdl.transforms.base.linking.expression.values

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class EngineFunctionEvaluatorsSpec extends AnyFlatSpec with Matchers{
  behavior of "read_lines"
  it should "return an empty array from an empty file" in {
    EngineFunctionEvaluators.split_lines("") shouldBe(List.empty)
  }
  it should "return a one element array with a one line file" in {
    EngineFunctionEvaluators.split_lines("one line") shouldBe(List("one line"))
  }
  it should "return a one element array with a one line file ending with a new line character" in {
    EngineFunctionEvaluators.split_lines("one line\n") shouldBe(List("one line"))
  }
  it should "return a one element array with a single new line character" in {
    EngineFunctionEvaluators.split_lines("\n") shouldBe(List(""))
  }
  it should "return a three element array with three new line characters" in {
    EngineFunctionEvaluators.split_lines("\n\n\n") shouldBe(List(""))
  }
  it should "return a three element array with three words and three double new line characters" in {
    EngineFunctionEvaluators.split_lines("a word\n\na word\n\none more word\n\n") shouldBe(List("a word", "", "a word", "", "one more word"))
  }
}
