package cromwell.core.simpleton

import cromwell.core.CallOutputs
import cromwell.core.simpleton.WomValueSimpleton._
import wom.graph.GraphNodePort.OutputPort
import wom.types._
import wom.values._

import scala.language.postfixOps


/**
  * Builds arbitrary `WomValues` from `WomValueSimpletons`.
  **/
object WomValueBuilder {

  /**
    * Looks for a WOM identifier possibly followed by a metacharacter and more stuff. A 'WOM identifier' is kind of a
    * made up term that encompasses anything preceding an open square brace or colon (used for encoding simpleton paths).
    *
    * Capture groups:
    *
    * <ol>
    * <li>A WOM identifier (nongreedy match of one or more of any character)</li>
    * <li>Either:
    *   <ol>
    *     <li>An open square bracket or colon followed by anything.</li>
    *     <li>Nothing at all, forcing the initial nongreedy match to consume everything.</li>
    *   </ol>
    * </li>
    * </ol>
    *
    */
  private val IdentifierAndPathPattern = raw"^(.+?)([\[:].*|)".r

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
    * metacharacter and more stuff.  Any metacharacters in map keys have been escaped by `WomValueSimpleton#escapeMeta`,
    * which is taken account by the regular expression.
    *
    * Capture groups:
    *
    * <ol>
    * <li>The map key, possibly including `WomValueSimpleton#escapeMeta` escaped metacharacters.  Does <b>not</b> include
    *     the leading colon.</li>
    * <li>Possibly a metacharacter and more trailing stuff</li>
    * </ol>
    */
  // Within the noncapturing `?:` group, this looks for an escaped metacharacter OR a non-metacharacter.
  private val MapElementPattern = raw"^:((?:\\[]\[:]|[^]\[:])+)(.*)".r

  // Group tuples by key using a Map with key type `K`.
  private def group[K](tuples: Traversable[(K, SimpletonComponent)]): Map[K, Traversable[SimpletonComponent]] = {
    tuples groupBy { case (i, _) => i } map { case (k, v) => k -> (v map { case (_, s) => s}) }
  }

  // Returns a tuple of the index into the outermost array and a `SimpletonComponent` whose path reflects the "descent"
  // into the array.  e.g. for a component
  // SimpletonComponent("[0][1]", v) this would return (0 -> SimpletonComponent("[1]", v)).
  private def descendIntoArray(component: SimpletonComponent): (Int, SimpletonComponent) = {
    component.path match { case ArrayElementPattern(index, more) => index.toInt -> component.copy(path = more)}
  }

  // Returns a tuple of the key into the outermost map and a `SimpletonComponent` whose path reflects the "descent"
  // into the map.  e.g. for a component
  // SimpletonComponent(":bar:baz", v) this would return ("bar" -> SimpletonComponent(":baz", v)).
  // Map keys are treated as Strings by this method, the caller must ultimately do the appropriate coercion to the
  // actual map key type.
  private def descendIntoMap(component: SimpletonComponent): (String, SimpletonComponent) = {
    component.path match { case MapElementPattern(key, more) => key.unescapeMeta -> component.copy(path = more)}
  }

  private implicit class EnhancedSimpletonComponents(val components: Traversable[SimpletonComponent]) extends AnyVal {
    def asArray: List[Traversable[SimpletonComponent]] = group(components map descendIntoArray).toList.sortBy(_._1).map(_._2)
    def asMap: Map[String, Traversable[SimpletonComponent]] = group(components map descendIntoMap)
    def asPrimitive: WomValue = components.head.value
    def asString: String = asPrimitive.valueString
  }

  private def toWomValue(outputType: WomType, components: Traversable[SimpletonComponent]): WomValue = {



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
    
    def toWomFile(components: Traversable[SimpletonComponent]) = {
      // If there's just one simpleton, it's a primitive (file or directory)
      if (components.size == 1) components.asPrimitive
      else {
        // Otherwise make a map of the components and detect the type of file from the class field
        val groupedListing = components.asMap
        
        def isClass(className: String) = {
          groupedListing.get(ClassKey)
          /* If the class field is in an array it will be prefixed with a ':', so check for that as well.
           * e.g: secondaryFiles[0]:class -> "File" 
           *      secondaryFiles[0]:value -> "file/path" 
           * would produce a Map(
           *  ":class" -> List(Simpleton("File")),
           *  ":value" -> List(Simpleton("file/path"))
           * )
           */
          .orElse(groupedListing.get(s":$ClassKey"))
          .map(_.asPrimitive.valueString)
          .contains(className)
        } 
          
        def isDirectory = isClass(WomValueSimpleton.DirectoryClass)
        def isFile = isClass(WomValueSimpleton.FileClass)

        if (isDirectory) toWomValue(WomMaybeListedDirectoryType, components)
        else if (isFile) toWomValue(WomMaybePopulatedFileType, components)
        else throw new IllegalArgumentException(s"There is no WomFile that can be built from simpletons: ${groupedListing.toList.mkString(", ")}")
      }
    }
    
    outputType match {
      case _: WomPrimitiveType =>
        components.asPrimitive
      case opt: WomOptionalType =>
        if (components.isEmpty) {
          WomOptionalValue(opt.memberType, None)
        } else {
          WomOptionalValue(toWomValue(opt.memberType, components))
        }
      case arrayType: WomArrayType =>
        WomArray(arrayType, components.asArray map { toWomValue(arrayType.memberType, _) })
      case mapType: WomMapType =>
        // map keys are guaranteed by WOM to be primitives, so the "coerceRawValue(..).get" is safe.
        WomMap(mapType, components.asMap map { case (k, ss) => mapType.keyType.coerceRawValue(k).get -> toWomValue(mapType.valueType, ss) })
      case pairType: WomPairType =>
        val groupedByLeftOrRight: Map[PairLeftOrRight, Traversable[SimpletonComponent]] = group(components map descendIntoPair)
        WomPair(toWomValue(pairType.leftType, groupedByLeftOrRight(PairLeft)), toWomValue(pairType.rightType, groupedByLeftOrRight(PairRight)))
      case WomObjectType =>
        // map keys are guaranteed by WOM to be primitives, so the "coerceRawValue(..).get" is safe.
        val map: Map[String, WomValue] = components.asMap map { case (k, ss) => k -> toWomValue(WomAnyType, ss) }
        WomObject(map)
      case composite: WomCompositeType =>
        val map: Map[String, WomValue] = components.asMap map { case (k, ss) => 
          val valueType = composite
            .typeMap
            .getOrElse(k, throw new RuntimeException(s"Field $k is not a declared field of composite type $composite. Cannot build a WomValue from the simpletons."))
          k -> toWomValue(valueType, ss) 
        }
        WomObject.withTypeUnsafe(map, composite)
      case WomMaybeListedDirectoryType =>
        val directoryValues = components.asMap

        val value = directoryValues.get("value").map(_.asString)
        val listing = directoryValues.get("listing")
          .map({ _.asArray.map(toWomFile).collect({ case womFile: WomFile => womFile }) })

        WomMaybeListedDirectory(value, listing)
      case WomMaybePopulatedFileType =>
        val populatedValues = components.asMap

        val value = populatedValues.get("value").map(_.asString)
        val checksum = populatedValues.get("checksum").map(_.asString)
        val size = populatedValues.get("size").map(_.asString.toLong)
        val format = populatedValues.get("format").map(_.asString)
        val contents = populatedValues.get("contents").map(_.asString)
        val secondaryFiles = populatedValues.get("secondaryFiles").toList.flatMap({ 
          _.asArray.map(toWomFile).collect({ case womFile: WomFile => womFile })
        })

        WomMaybePopulatedFile(
          valueOption = value,
          checksumOption = checksum,
          sizeOption = size,
          formatOption = format,
          contentsOption = contents,
          secondaryFiles = secondaryFiles
        )
      case coproductType: WomCoproductType =>
        // We don't currently record the actual type of the coproduct value so use the same heuristics as for Any.
        WomCoproductValue(coproductType, toWomValue(WomAnyType, components))

      case WomAnyType =>
        // Ok, we're going to have to guess, but the keys should give us some clues:
        if (components forall { component => MapElementPattern.findFirstMatchIn(component.path).isDefined }) {
          // Looks like an object
          toWomValue(WomObjectType, components)
        } else if (components forall { component => ArrayElementPattern.findFirstMatchIn(component.path).isDefined }) {
          // Looks like an array:
          toWomValue(WomArrayType(WomAnyType), components) match {
            case WomArray(_, values) => WomArray(values)
          }

        } else {
          // Treat it as a primitive:
          components collectFirst { case SimpletonComponent(_, v) => v } get
        }
    }
  }

  /** The two fields of the `SimpletonComponent` product type have types identical to those of the `WomValueSimpleton`
    * product type, but a `SimpletonComponent` is conceptually different from a `WomValueSimpleton`.
    * `SimpletonComponent` removes the outermost "name" piece of the `WomValueSimpleton` key, so `path` contains only the
    * path to the element.  e.g. for a `WomValueSimpleton` of
    *
    * {{{
    * WomValueSimpleton("foo:bar[0]", WomString("baz"))
    * }}}
    *
    * the corresponding `SimpletonComponent` would be
    *
    * {{{
    * SimpletonComponent(":bar[0]", WomString("baz"))
    * }}}
    */
  private case class SimpletonComponent(path: String, value: WomValue)

  def toJobOutputs(taskOutputs: Traversable[OutputPort], simpletons: Traversable[WomValueSimpleton]): CallOutputs = {
    CallOutputs(toWomValues(taskOutputs, simpletons))
  }

  def toWomValues(taskOutputs: Traversable[OutputPort], simpletons: Traversable[WomValueSimpleton]): Map[OutputPort, WomValue] = {

    def simpletonToComponent(name: String)(simpleton: WomValueSimpleton): SimpletonComponent = {
      SimpletonComponent(simpleton.simpletonKey.drop(name.length), simpleton.simpletonValue)
    }

    // This is meant to "rehydrate" simpletonized WomValues back to WomValues.  It is assumed that these WomValues were
    // "dehydrated" to WomValueSimpletons correctly. This code is not robust to corrupt input whatsoever.
    val types = taskOutputs map { o => o -> o.womType } toMap
    val simpletonsByOutputName = simpletons groupBy { _.simpletonKey match { case IdentifierAndPathPattern(i, _) => i } }
    val simpletonComponentsByOutputName: Map[String, Traversable[SimpletonComponent]] =
      simpletonsByOutputName map { case (name, ss) => name -> (ss map simpletonToComponent(name)) }
    types map { case (outputPort, outputType) => outputPort -> toWomValue(outputType, simpletonComponentsByOutputName.getOrElse(outputPort.internalName, Seq.empty))}
  }
}

