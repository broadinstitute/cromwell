package cwl

import common.validation.Validation._
import eu.timepit.refined.refineMV
import cwl.CwlDecoder.decodeCwlFile
import cwl.TestSetup.rootPath
import eu.timepit.refined.numeric.Positive
import org.scalatest.{FlatSpec, Matchers}
import wom.callable.Callable.InputDefinition
import wom.callable.{CallableTaskDefinition, RuntimeEnvironment}
import wom.expression.NoIoFunctionSet
import wom.graph.OptionalGraphInputNodeWithDefault
import wom.values.WomValue

class DirectorySpec extends FlatSpec with Matchers {

  behavior of "Directory"

  it should "dir_example" in {
    val cwl = decodeCwlFile(rootPath / "dir_example.cwl").value.unsafeRunSync.right.get
    val executable = cwl.womExecutable(AcceptAllRequirements, None, NoIoFunctionSet, strictValidation = false).right.get
    val call = executable.graph.calls.head
    val runtimeEnvironment = RuntimeEnvironment("output/path", "temp/path",refineMV[Positive](1), 2e10, 100, 100)
    val defaultCallInputs = executable.graph.nodes.collect({
      case oginwd: OptionalGraphInputNodeWithDefault =>
        val key: InputDefinition = call.inputDefinitionMappings.toMap.keys.find(
          _.localName == oginwd.identifier.localName
        ).get
        val value: WomValue = oginwd.default.evaluateValue(Map.empty, NoIoFunctionSet).toTry.get
        key -> value
    }).toMap
    val command = call.callable.asInstanceOf[CallableTaskDefinition].instantiateCommand(
      defaultCallInputs, NoIoFunctionSet, identity, runtimeEnvironment
    ).toTry.get
    command.commandString should be("""   'echo' 'exampledir'   """.trim)
  }

}
