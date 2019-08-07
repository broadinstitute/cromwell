package wdl

import better.files._
import org.scalatest.{FlatSpec, Matchers}
import spray.json._
import wdl.draft2.Draft2ResolvedImportBundle
import wdl.draft2.model._
import wom.ResolvedImportRecord

import scala.collection.immutable.ListMap

class WdlWiringSpec extends FlatSpec with Matchers {
  val testCases = File("src/test/cases")
  testCases.createDirectories()

  testCases.list.toSeq.filter(_.isDirectory) foreach { testDir =>
    val wdlFile = testDir / "test.wdl"
    if (!wdlFile.exists) fail(s"Expecting a 'test.wdl' file in directory 'cases/${testDir.name}'")
    def resolvers: Seq[Draft2ImportResolver] =
      Seq((relPath: String) => Draft2ResolvedImportBundle((testDir / relPath).contentAsString, ResolvedImportRecord((testDir / relPath).pathAsString)))
    val namespace = WdlNamespaceWithWorkflow.load(File(wdlFile.path).contentAsString, resolvers).get
    val wdlFileRelPath = File(".").relativize(wdlFile)

    expectedFullyQualifiedNames(testDir, namespace) foreach { case (fqn, expectedType) =>
      it should s"resolve FQN $fqn to object of type $expectedType in WDL file $wdlFileRelPath" in {
        val resolution = namespace.resolve(fqn)
        resolution.map(_.getClass.getSimpleName) shouldEqual Option(expectedType)
        resolution.map(_.fullyQualifiedName) shouldEqual Option(fqn)
      }
    }

    expectedInputs(testDir, namespace) foreach { case (fqn, womType) =>
      it should s"have $fqn (of type $womType) as an input in WDL file $wdlFileRelPath" in {
        val input = namespace.workflow.inputs.get(fqn)
        input should not be None
        input.map(_.womType.stableName) shouldEqual Option(womType)
      }
    }

    expectedFullyQualifiedNamesWithIndexScopes(testDir, namespace) foreach { case (fqn, expectedType) =>
      it should s"resolve FQN (with index scopes) $fqn to object of type $expectedType in WDL file $wdlFileRelPath" in {
        val resolution = namespace.resolve(fqn)
        resolution.map(_.getClass.getSimpleName) shouldEqual Option(expectedType)
        resolution.map(_.fullyQualifiedNameWithIndexScopes) shouldEqual Option(fqn)
      }
    }

    expectedParents(testDir, namespace) foreach { case (nodeFqn, parentFqn) =>
      it should s"compute parent of $nodeFqn to be $parentFqn in WDL file $wdlFileRelPath" in {
        val nodeResolution = namespace.resolve(nodeFqn)
        val parentResolution = parentFqn.flatMap(namespace.resolve)
        nodeResolution.flatMap(_.parent) shouldEqual parentResolution
      }
    }

    expectedChildren(testDir, namespace) foreach { case (nodeFqn, children) =>
      it should s"compute children of $nodeFqn to be ${children.map(_.fullyQualifiedName).mkString(", ")} in WDL file $wdlFileRelPath" in {
        val nodeResolution = namespace.resolve(nodeFqn)
        nodeResolution.map(_.children) shouldEqual Option(children)
      }
    }

    expectedUpstream(testDir, namespace) foreach { case (node, expectedUpstreamNodes) =>
      it should s"compute upstream nodes for FQN ${node.fullyQualifiedName} in WDL file $wdlFileRelPath" in {
        node.upstream shouldEqual expectedUpstreamNodes
      }
    }

    expectedUpstreamAncestry(testDir, namespace) foreach { case (node, expectedUpstreamAncestorNodes) =>
      it should s"compute full upstream ancestry for FQN ${node.fullyQualifiedName} in WDL file $wdlFileRelPath" in {
        node.upstreamAncestry shouldEqual expectedUpstreamAncestorNodes
      }
    }

    expectedDownstream(testDir, namespace) foreach { case (node, expectedDownstreamNodes) =>
      it should s"compute downstream nodes for FQN ${node.fullyQualifiedName} in WDL file $wdlFileRelPath" in {
        node.downstream shouldEqual expectedDownstreamNodes
      }
    }

    expectedAncestry(testDir, namespace) foreach { case (node, expectedAncestry) =>
      it should s"compute ancestry for FQN ${node.fullyQualifiedName} in WDL file $wdlFileRelPath" in {
        node.ancestry shouldEqual expectedAncestry
      }
    }
  }

  private def expectedInputs(testDir: File, namespace: WdlNamespaceWithWorkflow): Map[FullyQualifiedName, String] = {
    val expectedWorkflowInputsFile = testDir / "inputs.expectations"

    if (!expectedWorkflowInputsFile.exists) {
      val workflowInputs = namespace.workflow.inputs map { case (fqn, input) =>
        fqn -> JsString(input.womType.stableName)
      }
      val jsObject = JsObject(ListMap(workflowInputs.toSeq.sortBy(_._1): _*))
      expectedWorkflowInputsFile.write(jsObject.prettyPrint + "\n")
    }

    expectedWorkflowInputsFile.contentAsString.parseJson.asInstanceOf[JsObject].fields.asInstanceOf[Map[String, JsString]] map {
      case (k, v) => k -> v.value
    }
  }

  private def expectedParents(testDir: File, namespace: WdlNamespaceWithWorkflow): Map[FullyQualifiedName, Option[FullyQualifiedName]] = {
    val expectedParentsFile = testDir / "parents.expectations"

    if (!expectedParentsFile.exists) {
      val fqnsAndParent = namespace.descendants map { scope =>
        scope.fullyQualifiedName -> scope.parent.map(_.fullyQualifiedName).map(JsString(_)).getOrElse(JsNull)
      }
      val jsObject = JsObject(ListMap(fqnsAndParent.toSeq.sortBy(_._1): _*))
      expectedParentsFile.write(jsObject.prettyPrint + "\n")
    }

    expectedParentsFile.contentAsString.parseJson.asInstanceOf[JsObject].fields.asInstanceOf[Map[String, JsValue]] map {
      case (k, v: JsString) => k -> Option(v.value)
      case (k, _) => k -> None
    }
  }

  private def expectedChildren(testDir: File, namespace: WdlNamespaceWithWorkflow): Map[FullyQualifiedName, Seq[Scope]] = {
    val expectedChildrenFile = testDir / "children.expectations"

    if (!expectedChildrenFile.exists) {
      val fqnsAndChildren = namespace.descendants map { scope =>
        scope.fullyQualifiedName -> JsArray(scope.children.map(_.fullyQualifiedName).map(JsString(_)).toVector)
      }
      val jsObject = JsObject(ListMap(fqnsAndChildren.toSeq.sortBy(_._1): _*))
      expectedChildrenFile.write(jsObject.prettyPrint + "\n")
    }

    expectedChildrenFile.contentAsString.parseJson.asInstanceOf[JsObject].fields.asInstanceOf[Map[String, JsArray]] map {
      case (k, v) =>
        val children = v.elements.collect({ case s: JsString => s }).map(s => namespace.resolve(s.value).get)
        k -> children
    }
  }

  private def expectedFullyQualifiedNames(testDir: File, namespace: WdlNamespaceWithWorkflow): Map[FullyQualifiedName, String] = {
    val expectedFqnsAndClassFile = testDir / "fqn.expectations"

    if (!expectedFqnsAndClassFile.exists) {
      val fqnsAndClassType = namespace.descendants map { scope =>
        scope.fullyQualifiedName -> JsString(scope.getClass.getSimpleName)
      }
      val jsObject = JsObject(ListMap(fqnsAndClassType.toSeq.sortBy(_._1): _*))
      expectedFqnsAndClassFile.write(jsObject.prettyPrint + "\n")
    }

    expectedFqnsAndClassFile.contentAsString.parseJson.asInstanceOf[JsObject].fields.asInstanceOf[Map[String, JsString]] map {
      case (k, v) => k -> v.value
    }
  }

  private def expectedFullyQualifiedNamesWithIndexScopes(testDir: File, namespace: WdlNamespaceWithWorkflow): Map[FullyQualifiedName, String] = {
    val expectedFqnsAndClassFile = testDir / "fqn_index_scopes.expectations"

    if (!expectedFqnsAndClassFile.exists) {
      val fqnsAndClassType = namespace.descendants map { scope =>
        scope.fullyQualifiedNameWithIndexScopes -> JsString(scope.getClass.getSimpleName)
      }
      val jsObject = JsObject(ListMap(fqnsAndClassType.toSeq.sortBy(_._1): _*))
      expectedFqnsAndClassFile.write(jsObject.prettyPrint + "\n")
    }

    expectedFqnsAndClassFile.contentAsString.parseJson.asInstanceOf[JsObject].fields.asInstanceOf[Map[String, JsString]] map {
      case (k, v) => k -> v.value
    }
  }

  private def expectedAncestry(testDir: File, namespace: WdlNamespaceWithWorkflow): Map[Scope, Seq[Scope]] = {
    val expectedAncestryFile = testDir / "ancestry.expectations"

    if (!expectedAncestryFile.exists) {
      val ancestryFqns = namespace.descendants map { scope =>
        scope.fullyQualifiedName -> JsArray(scope.ancestry.toVector.map(_.fullyQualifiedName).map(JsString(_)))
      }
      val jsObject = JsObject(ListMap(ancestryFqns.toSeq.sortBy(_._1): _*))
      expectedAncestryFile.write(jsObject.prettyPrint + "\n")
    }

    expectedAncestryFile.contentAsString.parseJson.asInstanceOf[JsObject].fields.asInstanceOf[Map[String, JsArray]] map {
      case (k, v) =>
        val expectedAncestry = v.elements.asInstanceOf[Vector[JsString]].map(n => namespace.resolve(n.value).get)
        val resolvedFqn = namespace.resolve(k).get
        resolvedFqn -> expectedAncestry
    }
  }

  private def expectedUpstreamAncestry(testDir: File, namespace: WdlNamespaceWithWorkflow): Map[WdlGraphNode, Set[Scope]] = {
    val expectedUpstreamFile = testDir / "upstreamAncestry.expectations"

    if (!expectedUpstreamFile.exists) {
      val upstreamFqns = namespace.descendants.collect({ case n: WdlGraphNode => n }) map { node =>
        node.fullyQualifiedName -> JsArray(node.upstreamAncestry.toVector.map(_.fullyQualifiedName).sorted.map(JsString(_)))
      }
      val jsObject = JsObject(ListMap(upstreamFqns.toSeq.sortBy(_._1): _*))
      expectedUpstreamFile.write(jsObject.prettyPrint + "\n")
    }

    expectedUpstreamFile.contentAsString.parseJson.asInstanceOf[JsObject].fields.asInstanceOf[Map[String, JsArray]] map {
      case (k, v) =>
        val expectedUpstreamAncestors = v.elements.asInstanceOf[Vector[JsString]].map(n => namespace.resolve(n.value).get).toSet
        val resolvedFqn = namespace.resolve(k).get.asInstanceOf[WdlGraphNode]
        resolvedFqn -> expectedUpstreamAncestors
    }
  }

  private def expectedUpstream(testDir: File, namespace: WdlNamespaceWithWorkflow): Map[WdlGraphNode, Set[Scope]] = {
    val expectedUpstreamFile = testDir / "upstream.expectations"

    if (!expectedUpstreamFile.exists) {
      val upstreamFqns = namespace.descendants.collect({ case n: WdlGraphNode => n }) map { node =>
        node.fullyQualifiedName -> JsArray(node.upstream.toVector.map(_.fullyQualifiedName).sorted.map(JsString(_)))
      }
      val jsObject = JsObject(ListMap(upstreamFqns.toSeq.sortBy(_._1): _*))
      expectedUpstreamFile.write(jsObject.prettyPrint + "\n")
    }

    expectedUpstreamFile.contentAsString.parseJson.asInstanceOf[JsObject].fields.asInstanceOf[Map[String, JsArray]] map {
      case (k, v) =>
        val expectedUpstream = v.elements.asInstanceOf[Vector[JsString]].map(n => namespace.resolve(n.value).get).toSet
        val resolvedFqn = namespace.resolve(k).get.asInstanceOf[WdlGraphNode]
        resolvedFqn -> expectedUpstream
    }
  }

  private def expectedDownstream(testDir: File, namespace: WdlNamespaceWithWorkflow): Map[WdlGraphNode, Set[Scope]] = {
    val expectedDownstreamFile = testDir / "downstream.expectations"

    if (!expectedDownstreamFile.exists) {
      val downstreamFqns = namespace.descendants.collect({ case n: WdlGraphNode => n }) map { node =>
        node.fullyQualifiedName -> JsArray(node.downstream.toVector.map(_.fullyQualifiedName).sorted.map(JsString(_)))
      }
      val jsObject = JsObject(ListMap(downstreamFqns.toSeq.sortBy(_._1): _*))
      expectedDownstreamFile.write(jsObject.prettyPrint + "\n")
    }

    expectedDownstreamFile.contentAsString.parseJson.asInstanceOf[JsObject].fields.asInstanceOf[Map[String, JsArray]] map { case (k, v) =>
      val expectedDownstream = v.elements.asInstanceOf[Vector[JsString]].map(n => namespace.resolve(n.value).get).toSet
      val resolvedFqn = namespace.resolve(k).get.asInstanceOf[WdlGraphNode]
      resolvedFqn -> expectedDownstream
    }
  }
}
