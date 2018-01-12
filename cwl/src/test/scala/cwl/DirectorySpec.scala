package cwl

import common.validation.Validation._
import cwl.CwlDecoder.decodeAllCwl
import cwl.TestSetup.rootPath
import org.scalatest.{FlatSpec, Matchers}
import wom.callable.Callable.InputDefinition
import wom.callable.{CallableTaskDefinition, RuntimeEnvironment}
import wom.expression.PlaceholderIoFunctionSet
import wom.graph.OptionalGraphInputNodeWithDefault
import wom.values.WomValue

class DirectorySpec extends FlatSpec with Matchers {

  behavior of "Directory"

  it should "dir_example" in {
    val cwl = decodeAllCwl(rootPath / "dir_example.cwl").value.unsafeRunSync.right.get
    val executable = cwl.womExecutable(AcceptAllRequirements, None).right.get
    val call = executable.graph.calls.head
    val runtimeEnvironment = RuntimeEnvironment("output/path", "temp/path", 1, 2e10, 100, 100)
    val defaultCallInputs = executable.graph.nodes.collect({
      case oginwd: OptionalGraphInputNodeWithDefault =>
        val key: InputDefinition = call.inputDefinitionMappings.toMap.keys.find(
          _.localName == oginwd.identifier.localName
        ).get
        val value: WomValue = oginwd.default.evaluateValue(Map.empty, PlaceholderIoFunctionSet).toTry.get
        key -> value
    }).toMap
    val command = call.callable.asInstanceOf[CallableTaskDefinition].instantiateCommand(
      defaultCallInputs, PlaceholderIoFunctionSet, identity, runtimeEnvironment
    ).toTry.get
    command.commandString should be("""   "echo" exampledir   """.trim)
  }

}
