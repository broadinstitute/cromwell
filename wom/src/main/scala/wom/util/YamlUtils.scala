package wom.util

import java.io.StringReader

import com.typesafe.config.ConfigException.BadValue
import com.typesafe.config.{Config, ConfigFactory}
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.NonNegative
import eu.timepit.refined.refineV
import io.circe.{Json, ParsingFailure}
import net.ceedubs.ficus.Ficus._
import net.ceedubs.ficus.readers.ValueReader

import scala.collection.JavaConverters._

object YamlUtils {

  private[util] implicit val refinedNonNegativeReader: ValueReader[Int Refined NonNegative] = {
    (config: Config, path: String) => {
      val int = config.getInt(path)
      refineV[NonNegative](int) match {
        case Left(error) => throw new BadValue(path, error)
        case Right(refinedInt) => refinedInt
      }
    }
  }

  private val defaultMaxNodes = ConfigFactory.load().as[Int Refined NonNegative]("yaml.max-nodes")

  /**
    * Parses the yaml, detecting loops, and a maximum number of nodes.
    *
    * See: https://en.wikipedia.org/wiki/Billion_laughs_attack
    *
    * @param yaml     The yaml string.
    * @param maxNodes Maximum number of yaml nodes to parse.
    * @return The parsed yaml.
    */
  def parse(yaml: String, maxNodes: Int Refined NonNegative = defaultMaxNodes): Either[ParsingFailure, Json] = {
    try {
      val snakeYamlParser = new org.yaml.snakeyaml.Yaml()
      val parsed = snakeYamlParser.load[AnyRef](new StringReader(yaml))
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
