package cromwell.languages

import wom.executable.Executable
import wom.graph.GraphNodePort.OutputPort
import wom.values.WomValue

/**
  * @param executable The WOM executable produced
  * @param womValueInputs The inputs to the executable, as interpreted from the inputs file
  * @param importedFileContent For reference, a mapping of {URI => file content} from any imported files during the WOMification process.
  */
final case class ValidatedWomNamespace(executable: Executable,
                                       womValueInputs: Map[OutputPort, WomValue],
                                       importedFileContent: Map[String, String]
)
