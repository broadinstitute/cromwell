package wom.executable

import wom.graph.Graph
import wom.graph.GraphNodePort.OutputPort
import wom.values.WomValue

case class ValidatedWomNamespace(executable: Executable, graph: Graph, womValueInputs: Map[OutputPort, WomValue])
