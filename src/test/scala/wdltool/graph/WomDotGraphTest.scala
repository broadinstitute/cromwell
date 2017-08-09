package wdltool.graph

import java.util.concurrent.atomic.AtomicInteger
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wom.graph.Graph

trait WomDotGraphTest extends FlatSpec with Matchers {

  def cases: List[WomDotGraphTestCase]

  behavior of "womgraph"

  def tests() = {
    cases foreach { testCase =>
      it should s"draw the right DOT graph for ${testCase.name}" in {
        val womGraph = new WomGraph(testCase.name, testCase.graph)

        val standardResult = standardizifyResult(womGraph.digraphDot)
        val standardExpectation = standardizifyResult(testCase.dotExpectation)

        standardResult should be(standardExpectation)
      }
    }
  }
  
  def standardizifyResult(dot: String) = {
    val regex = "(NODE|PORT)-?[0-9]*".r
    val matches = regex.findAllIn(dot).toList.distinct.map(_.toString)
    val currentMatchId = new AtomicInteger(0)
    def standardize(m: String) = {
      m.substring(0, 4) + currentMatchId.getAndIncrement().toString
    }
    def foldFunction(acc: String, m: String): String = acc.replaceAll(m, standardize(m))

    matches.foldLeft(dot)(foldFunction).replaceAll("\n[\\s]*", "\n")

  }
}


final case class WomDotGraphTestCase(name: String, graph: Graph, dotExpectation: String)
