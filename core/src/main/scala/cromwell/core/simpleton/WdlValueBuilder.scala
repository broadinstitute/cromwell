package cromwell.core.simpleton

import wdl4s.TaskOutput
import wdl4s.types._
import wdl4s.values.{WdlArray, WdlMap, WdlOptionalValue, WdlPair, WdlValue}

import scala.language.postfixOps
import cromwell.core.{CallOutputs, JobOutput}
import cromwell.core.simpleton.WdlValueSimpleton._


/**
  * Builds arbitrary `WdlValues` from `WdlValueSimpletons`.
  **/
object WdlValueBuilder {

  /**
    * Looks for a WDL identifier possibly followed by a metacharacter and more stuff.
    *
    * Capture groups:
    *
    * <ol>
    * <li>A WDL identifier</li>
    * <li>Possibly a metacharacter and more stuff after the WDL identifier</li>
    * </ol>
    */
  private val IdentifierAndPathPattern = "^([a-zA-Z][a-zA-Z0-9_]*)(.*)".r

  /**
    * Looks for an array element reference: square braces surrounding digits, possibly followed by a metacharacter and
    * more stuff.
    *
    * Capture groups:
    *
    * <ol>
    * <li>The array subscript. Does <b>not</b> include the enclosing square braces.</li>
    * <li>Possibly a metacharacter and more stuff after the square brace enclosed digits</li>
    * </ol>
    */
  private val ArrayElementPattern = raw"^\[(\d+)\](.*)".r

  /**
    * Looks for a map element reference: a colon followed by one or more non-metacharacters, possibly followed by a
    * metacharacter and more stuff.  Any metacharacters in map keys have been escaped by `WdlValueSimpleton#escapeMeta`,
    * which is taken account by the regular expression.
    *
    * Capture groups:
    *
    * <ol>
    * <li>The map key, possibly including `WdlValueSimpleton#escapeMeta` escaped metacharacters.  Does <b>not</b> include
    *     the leading colon.</li>
    * <li>Possibly a metacharacter and more trailing stuff</li>
    * </ol>
    */
  // Within the noncapturing `?:` group, this looks for an escaped metacharacter OR a non-metacharacter.
  private val MapElementPattern = raw"^:((?:\\[]\[:]|[^]\[:])+)(.*)".r

  private def toWdlValue(outputType: WdlType, components: Traversable[SimpletonComponent]): WdlValue = {

    // Returns a tuple of the index into the outermost array and a `SimpletonComponent` whose path reflects the "descent"
    // into the array.  e.g. for a component
    // SimpletonComponent("[0][1]", v) this would return (0 -> SimpletonComponent("[1]", v)).
    def descendIntoArray(component: SimpletonComponent): (Int, SimpletonComponent) = {
      component.path match { case ArrayElementPattern(index, more) => index.toInt -> component.copy(path = more)}
    }

    // Returns a tuple of the key into the outermost map and a `SimpletonComponent` whose path reflects the "descent"
    // into the map.  e.g. for a component
    // SimpletonComponent(":bar:baz", v) this would return ("bar" -> SimpletonComponent(":baz", v)).
    // Map keys are treated as Strings by this method, the caller must ultimately do the appropriate coercion to the
    // actual map key type.
    def descendIntoMap(component: SimpletonComponent): (String, SimpletonComponent) = {
      component.path match { case MapElementPattern(key, more) => key.unescapeMeta -> component.copy(path = more)}
    }

    // Returns a tuple of the key into the pair (i.e. left or right) and a `SimpletonComponent` whose path reflects the "descent"
    // into the pair.  e.g. for a component
    // SimpletonComponent(":left:foo", someValue) this would return (PairLeft -> SimpletonComponent(":baz", someValue)).
    sealed trait PairLeftOrRight
    case object PairLeft extends PairLeftOrRight
    case object PairRight extends PairLeftOrRight
    def descendIntoPair(component: SimpletonComponent): (PairLeftOrRight, SimpletonComponent) = {
      component.path match {
        case MapElementPattern("left", more) => PairLeft -> component.copy(path = more)
        case MapElementPattern("right", more) => PairRight -> component.copy(path = more)
      }
    }

    // Group tuples by key using a Map with key type `K`.
    def group[K](tuples: Traversable[(K, SimpletonComponent)]): Map[K, Traversable[SimpletonComponent]] = {
      tuples groupBy { case (i, _) => i } mapValues { _ map { case (i, s) => s} }
    }

    outputType match {
      case _: WdlPrimitiveType => components collectFirst { case SimpletonComponent(_, v) => v } get
      case opt: WdlOptionalType =>
        if (components.isEmpty) {
          WdlOptionalValue(opt.memberType, None)
        } else {
          WdlOptionalValue(toWdlValue(opt.memberType, components))
        }
      case arrayType: WdlArrayType =>
        val groupedByArrayIndex: Map[Int, Traversable[SimpletonComponent]] = group(components map descendIntoArray)
        WdlArray(arrayType, groupedByArrayIndex.toList.sortBy(_._1) map { case (_, s) => toWdlValue(arrayType.memberType, s) })
      case mapType: WdlMapType =>
        val groupedByMapKey: Map[String, Traversable[SimpletonComponent]] = group(components map descendIntoMap)
        // map keys are guaranteed by the WDL spec to be primitives, so the "coerceRawValue(..).get" is safe.
        WdlMap(mapType, groupedByMapKey map { case (k, ss) => mapType.keyType.coerceRawValue(k).get -> toWdlValue(mapType.valueType, ss) })
      case pairType: WdlPairType =>
        val groupedByLeftOrRight: Map[PairLeftOrRight, Traversable[SimpletonComponent]] = group(components map descendIntoPair)
        WdlPair(toWdlValue(pairType.leftType, groupedByLeftOrRight(PairLeft)), toWdlValue(pairType.rightType, groupedByLeftOrRight(PairRight)))
    }
  }

  /** The two fields of the `SimpletonComponent` product type have types identical to those of the `WdlValueSimpleton`
    * product type, but a `SimpletonComponent` is conceptually different from a `WdlValueSimpleton`.
    * `SimpletonComponent` removes the outermost "name" piece of the `WdlValueSimpleton` key, so `path` contains only the
    * path to the element.  e.g. for a `WdlValueSimpleton` of
    *
    * {{{
    * WdlValueSimpleton("foo:bar[0]", WdlString("baz"))
    * }}}
    *
    * the corresponding `SimpletonComponent` would be
    *
    * {{{
    * SimpletonComponent(":bar[0]", WdlString("baz"))
    * }}}
    */
  private case class SimpletonComponent(path: String, value: WdlValue)

  def toJobOutputs(taskOutputs: Traversable[TaskOutput], simpletons: Traversable[WdlValueSimpleton]): CallOutputs = {
    toWdlValues(taskOutputs, simpletons) mapValues JobOutput.apply
  }

  def toWdlValues(taskOutputs: Traversable[TaskOutput], simpletons: Traversable[WdlValueSimpleton]): Map[String, WdlValue] = {

    def simpletonToComponent(name: String)(simpleton: WdlValueSimpleton): SimpletonComponent = {
      SimpletonComponent(simpleton.simpletonKey.drop(name.length), simpleton.simpletonValue)
    }

    // This is meant to "rehydrate" simpletonized WdlValues back to WdlValues.  It is assumed that these WdlValues were
    // "dehydrated" to WdlValueSimpletons correctly. This code is not robust to corrupt input whatsoever.
    val types = taskOutputs map { o => o.unqualifiedName -> o.wdlType } toMap
    val simpletonsByOutputName = simpletons groupBy { _.simpletonKey match { case IdentifierAndPathPattern(i, _) => i } }
    val simpletonComponentsByOutputName = simpletonsByOutputName map { case (name, ss) => name -> (ss map simpletonToComponent(name)) }
    types map { case (name, outputType) => name -> toWdlValue(outputType, simpletonComponentsByOutputName(name))}
  }
}

