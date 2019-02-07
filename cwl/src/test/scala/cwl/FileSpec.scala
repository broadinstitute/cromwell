package cwl

import common.validation.Validation._
import eu.timepit.refined.refineMV
import cwl.CwlDecoder.decodeCwlFile
import cwl.TestSetup.rootPath
import eu.timepit.refined.numeric.Positive
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import wom.callable.Callable.InputDefinition
import wom.callable.{CallableTaskDefinition, RuntimeEnvironment}
import wom.graph.OptionalGraphInputNodeWithDefault
import wom.values.WomValue

class FileSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "File"

  private val fileTests = Table(
    ("description", "filePath", "ioFunctionSet", "expectedCommand"),
    ("get the basename of a file", "file_example.cwl", LocalIoFunctionSet, """   'echo' 'example.txt'   """.trim),
    ("get the size of a file", "file_size_example.cwl", LocalIoFunctionSet, """   'echo' '12'   """.trim)
  )

  forAll(fileTests) { (description, filePath, ioFunctionSet, expectedCommand) =>
    it should description in {
      val cwl = decodeCwlFile(rootPath / filePath).value.unsafeRunSync.right.get
      val executable = cwl.womExecutable(AcceptAllRequirements, None, ioFunctionSet, strictValidation = false).right.get
      val call = executable.graph.calls.head
      val runtimeEnvironment = RuntimeEnvironment("output/path", "temp/path", refineMV[Positive](1), 2e10, 100, 100)
      val defaultCallInputs = executable.graph.nodes.collect({
        case oginwd: OptionalGraphInputNodeWithDefault =>
          val key: InputDefinition = call.inputDefinitionMappings.toMap.keys.find(
            _.localName == oginwd.identifier.localName
          ).get
          val value: WomValue = key
            .valueMapper(ioFunctionSet)(oginwd.default.evaluateValue(Map.empty, ioFunctionSet).toTry.get)
            .value.unsafeRunSync()
            .right.get

          key -> value
      }).toMap
      val commandEither = call.callable.asInstanceOf[CallableTaskDefinition].instantiateCommand(
        defaultCallInputs, ioFunctionSet, identity, runtimeEnvironment
      ).toEither
      val command = commandEither.right.get.commandString
      command should be(expectedCommand)
    }
  }

}
