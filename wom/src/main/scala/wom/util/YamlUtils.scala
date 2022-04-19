package wom.util

import java.io.StringReader
import java.util
import com.typesafe.config.ConfigException.BadValue
import com.typesafe.config.{Config, ConfigFactory}
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.NonNegative
import eu.timepit.refined.refineV
import io.circe.{Json, ParsingFailure}
import net.ceedubs.ficus.Ficus._
import net.ceedubs.ficus.readers.ValueReader
import org.yaml.snakeyaml.LoaderOptions
import org.yaml.snakeyaml.comments.CommentLine
import org.yaml.snakeyaml.composer.Composer
import org.yaml.snakeyaml.constructor.SafeConstructor
import org.yaml.snakeyaml.nodes.{MappingNode, Node, NodeTuple}
import org.yaml.snakeyaml.parser.ParserImpl
import org.yaml.snakeyaml.reader.StreamReader
import org.yaml.snakeyaml.resolver.Resolver

import scala.jdk.CollectionConverters._

object YamlUtils {

  /**
    * Parses the yaml, detecting loops, and a maximum number of nodes.
    *
    * See: https://en.wikipedia.org/wiki/Billion_laughs_attack
    *
    * @param yaml     The yaml string.
    * @param maxNodes Maximum number of yaml nodes to parse.
    * @param maxDepth Maximum nested depth of yaml to parse.
    * @return The parsed yaml.
    */
  def parse(yaml: String,
            maxNodes: Int Refined NonNegative = defaultMaxNodes,
            maxDepth: Int Refined NonNegative = defaultMaxDepth
           ): Either[ParsingFailure, Json] = {
    try {
      val yamlConstructor = new SafeConstructor()
      val yamlComposer = new MaxDepthComposer(yaml, maxDepth)
      yamlConstructor.setComposer(yamlComposer)
      val parsed = yamlConstructor.getSingleData(classOf[AnyRef])

      // Use identity comparisons to check if two nodes are the same and do NOT recurse into them checking equality
      val identityHashMap = new java.util.IdentityHashMap[AnyRef, java.lang.Boolean]()
      // Since we don't actually need the values, wrap the Map in a Set
      val identitySet = java.util.Collections.newSetFromMap(identityHashMap)
      searchForOversizedYaml(parsed, identitySet, maxNodes, new Counter)
      io.circe.yaml.parser.parse(yaml)
    } catch {
      case exception: Exception =>
        Left(ParsingFailure(exception.getMessage, exception))
    }
  }

  private[util] implicit val refinedNonNegativeReader: ValueReader[Int Refined NonNegative] = {
    (config: Config, path: String) => {
      val int = config.getInt(path)
      refineV[NonNegative](int) match {
        case Left(error) => throw new BadValue(path, error)
        case Right(refinedInt) => refinedInt
      }
    }
  }

  private val yamlConfig = ConfigFactory.load().getConfig("yaml")
  private val defaultMaxNodes = yamlConfig.as[Int Refined NonNegative]("max-nodes")
  private val defaultMaxDepth = yamlConfig.as[Int Refined NonNegative]("max-depth")

  // Effectively turning off alias limits, inspired by
  // https://github.com/jenkinsci/configuration-as-code-plugin/issues/1374
  // https://github.com/jenkinsci/configuration-as-code-plugin/pull/1375
  // Cromwell still error checks and produces a hopefully more useful message than the
  // org.yaml.snakeyaml.error.YAMLException: Number of aliases for non-scalar nodes exceeds the specified max=50
  // that comes from the SnakeYAML library.
  private val loaderOptions = new LoaderOptions()
  loaderOptions.setMaxAliasesForCollections(Integer.MAX_VALUE)

  /** Extends SnakeYaml's Composer checking for a maximum depth before a StackOverflowError occurs. */
  private class MaxDepthComposer(yaml: String, maxDepth: Int Refined NonNegative)
    extends Composer(
      new ParserImpl(new StreamReader(new StringReader(yaml))),
      new Resolver(),
      loaderOptions
    ) {

    private val depth = new Counter

    private def checkDepth[A](f: => A): A = {
      depth.count += 1
      if (depth.count > maxDepth.value)
        throw new IllegalArgumentException(s"Parsing halted at node depth $maxDepth")
      val result = f
      depth.count -= 1
      result
    }

    override def composeScalarNode(anchor: String, blockComments: util.List[CommentLine]): Node = {
      checkDepth(super.composeScalarNode(anchor, blockComments))
    }

    override def composeSequenceNode(anchor: String): Node = {
      checkDepth(super.composeSequenceNode(anchor))
    }

    override def composeMappingNode(anchor: String): Node = {
      checkDepth(super.composeMappingNode(anchor))
    }

    override def composeMappingChildren(children: util.List[NodeTuple], node: MappingNode): Unit = {
      checkDepth(super.composeMappingChildren(children, node))
    }

    override def composeKeyNode(node: MappingNode): Node = {
      checkDepth(super.composeKeyNode(node))
    }

    override def composeValueNode(node: MappingNode): Node = {
      checkDepth(super.composeValueNode(node))
    }
  }

  /** A "pointer" reference to a mutable count. */
  private class Counter {
    var count = 0L
  }

  /**
    * Looks for loops and large documents in yaml parsed by SnakeYaml.
    *
    * Possibly can be refactored if/when Circe switches to SnakeYaml-Engine (aka Yaml 1.2), as SY-E disables recursive
    * key by default:
    *
    * - https://github.com/circe/circe-yaml/issues/46
    * - https://bitbucket.org/asomov/snakeyaml/issues/432/aggressive-yaml-anchors-causing
    * - https://bitbucket.org/asomov/snakeyaml-engine/commits/7573a8b4d551fc84521f4ac1234a361bfbc96698
    *
    * @param node        Current node to be evaluated.
    * @param identitySet Previously seen nodes during this branch of cycle checking.
    * @param maxNodes    The maximum number of nodes allowed to be traversed.
    * @param counter     A counter tracking to the total number of nodes traversed during the entire traversal.
    */
  private def searchForOversizedYaml(node: AnyRef,
                                     identitySet: java.util.Set[AnyRef],
                                     maxNodes: Int Refined NonNegative,
                                     counter: Counter): Unit = {
    if (!identitySet.add(node)) {
      throw new IllegalArgumentException("Loop detected")
    }

    counter.count += 1
    if (counter.count > maxNodes.value) {
      throw new IllegalArgumentException(s"Loop detection halted at $maxNodes nodes")
    }

    node match {
      case iterable: java.lang.Iterable[AnyRef]@unchecked =>
        iterable.asScala foreach {
          searchForOversizedYaml(_, identitySet, maxNodes, counter)
        }
      case map: java.util.Map[AnyRef, AnyRef]@unchecked =>
        map.asScala foreach {
          case (key, value) =>
            searchForOversizedYaml(key, identitySet, maxNodes, counter)
            searchForOversizedYaml(value, identitySet, maxNodes, counter)
        }
      case _ => /* ignore scalars, only loop through Yaml sequences and mappings: https://yaml.org/spec/1.1/#id861435 */
    }

    identitySet.remove(node)
    ()
  }
}
