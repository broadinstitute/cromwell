package cwl

import cwl.CommandLineTool.{CommandBindingSortingKey, CommandPartsList, SortKeyAndCommandPart, StringOrInt}
import cwl.CwlType.CwlType
import cwl.MyriadInputTypeToSortedCommandParts.CommandPartBuilder
import cwl.command.ParentName
import shapeless.{Coproduct, Poly1}
import wom.values.{WomArray, WomObjectLike, WomOptionalValue, WomValue}

/**
  * Poly1 to fold over a MyriadInputType and create a List of SortingKeyAndCommandPart based on the rules described here:
  * http://www.commonwl.org/v1.0/CommandLineTool.html#Input_binding
  */
object MyriadInputTypeToSortedCommandParts extends Poly1 {
  import Case._

  /**
    * CommandLineBinding: binding of the element being processed.
    *   It's possible for an input not to have an "inputBinding", but to have a complex type which itself has inputBindings for its fields or items.
    *   e.g:
    *     - id: arrayInput
    *       type:
    *         type: array
    *         items: File
    *         inputBinding: { prefix: "-YYY" }
    *         
    * Here the array type has an inputBinding but the input itself (arrayInput) does not.
    * This would yield on the command line something like "-YYY arrayItem0 -YYY arrayItem1 ..."
    * 
    * On the other hand:
    *     - id: arrayInput
    *       type:
    *         type: array
    *         items: File
    *         inputBinding: { prefix: "-YYY" }
    *       inputBinding: { prefix: "-XXX" }
    *       
    * In this case the input itself has a binding as well, which would yield something like "-XXX -YYY arrayItem0 -YYY arrayItem1 ..."
    * 
    * For arrays specifically, the presence of an "itemSeparator" and / or "valueFrom" in the input binding will change the way the array is processed,
    * which is why we need to know about it in this Poly1 that folds over the input types.
    *   
    * WomValue: value bound to the element (in the cases above it would be the WomArray)
    * CommandBindingSortingKey: The current sorting key. Because we might need to recurse in nested types we need to propagate the key as we do.
    */
  type CommandPartBuilder = (Option[InputCommandLineBinding], WomValue, CommandBindingSortingKey, ExpressionLib) => CommandPartsList

  implicit def m: Aux[MyriadInputInnerType, CommandPartBuilder] = at[MyriadInputInnerType] { innerType => {
    case (binding, womValue, sortingKey, expressionLib) =>
      innerType.fold(MyriadInputInnerTypeToSortedCommandParts).apply(binding, womValue, sortingKey, expressionLib)
  }
  }

  implicit def am: Aux[Array[MyriadInputInnerType], CommandPartBuilder] = at[Array[MyriadInputInnerType]] { types => {
    case (binding, womValue, sortingKey, expressionLib) =>
      types.toList.flatMap(_.fold(MyriadInputInnerTypeToSortedCommandParts).apply(binding, womValue, sortingKey, expressionLib))
  }
  }
}

object MyriadInputInnerTypeToSortedCommandParts extends Poly1 {

  import Case._

  def ex(component: String) = throw new RuntimeException(s"input type $component not yet supported by WOM!")

  // Primitive type: we just need to create a command part from the binding if there's one here.
  implicit def ct: Aux[CwlType, CommandPartBuilder] = at[CwlType] { _ => { case (binding, womValue, key, expressionLib) => binding.toList.map(_.toCommandPart(key, womValue, expressionLib)) } }

  implicit def irs: Aux[InputRecordSchema, CommandPartBuilder] = at[InputRecordSchema] { irs => {
    def go: CommandPartBuilder = {
      case (inputBinding, WomOptionalValue(_, Some(value)), sortingKey, expressionLib) =>
        go(inputBinding, value, sortingKey, expressionLib)
      case (inputBinding, objectLike: WomObjectLike, sortingKey, expressionLib) =>
        // If there's an input binding, make a SortKeyAndCommandPart for it
        val fromInputBinding: Option[SortKeyAndCommandPart] = inputBinding.map(_.toCommandPart(sortingKey, objectLike, expressionLib))

        // Go over the fields an fold over their type
        irs.fields.toList.flatten
          .foldLeft(CommandPartsList.empty)({
            case (currentMap, field) =>
              // Parse the name to get a clean id
              val parsedName = FullyQualifiedName(field.name)(ParentName.empty).id

              // The field name needs to be added to the key after the input binding (as per the spec)
              // Also start from the key from the input binding if there was one
              val fieldSortingKey =
              fromInputBinding.map(_.sortingKey).getOrElse(sortingKey)
                .append(field.inputBinding, Coproduct[StringOrInt](parsedName))

              // Get the value of this field from the objectLike
              // TODO: Could this fail ?
              val innerValue = objectLike.values(parsedName)

              val fromType = field.`type`.fold(MyriadInputTypeToSortedCommandParts).apply(field.inputBinding, innerValue, fieldSortingKey.asNewKey, expressionLib)

              currentMap ++ fromType
          }) ++ fromInputBinding
      case (_, other, _, _) => ex(s"Value $other cannot be used for an input of type InputRecordSchema")
    }

    go
  }}

  implicit def ies: Aux[InputEnumSchema, CommandPartBuilder] = at[InputEnumSchema] { iesType =>
    throw new RuntimeException(s"input type $iesType not yet supported by WOM!")
  }

  implicit def ias: Aux[InputArraySchema, CommandPartBuilder] = at[InputArraySchema] { ias =>
    def go: CommandPartBuilder = {
      case (inputBinding, WomOptionalValue(_, Some(value)), sortingKey, expressionLib) =>
        go(inputBinding, value, sortingKey, expressionLib)
      case (inputBinding, WomArray.WomArrayLike(womArray), sortingKey, expressionLib) =>

        // If there's an input binding, make a SortKeyAndCommandPart for it
        val fromInputBinding: Option[SortKeyAndCommandPart] = inputBinding.map(_.toCommandPart(sortingKey, womArray, expressionLib))

        // Now depending on whether we have an itemSeparator and/or valueFrom or not, we're going to recurse over each element of the array (or not).
        // See http://www.commonwl.org/v1.0/CommandLineTool.html#CommandLineBinding
        if (inputBinding.flatMap(_.itemSeparator).isDefined || inputBinding.flatMap(_.valueFrom).isDefined) {
          // If there's an item separator or a valueFrom we can stop here.
          // When the command part is instantiated (see CommandLineBindingCommandPart) it will evaluate the valueFrom (if defined) and join the items together (if there's an itemSeparator).
          fromInputBinding.toList
        } else {
          // If neither valueFrom nor itemSeparator were defined, we need to process each item of the array
          womArray.value.zipWithIndex.foldLeft(CommandPartsList.empty)({
            case (currentMap, (item, index)) =>
              // Update the sorting key with the binding position (if any), add the index
              val itemSortingKey =
              // The index needs to be added to the key after the input binding (as per the spec)
              // Also start from the key from the input binding if there was one
                fromInputBinding.map(_.sortingKey).getOrElse(sortingKey)
                  .append(ias.inputBinding, Coproduct[StringOrInt](index))

              // Fold over the item type fo each array element
              val fromType = ias.items.fold(MyriadInputTypeToSortedCommandParts).apply(ias.inputBinding, item, itemSortingKey.asNewKey, expressionLib)
              currentMap ++ fromType
          }) ++ fromInputBinding
        }

      case (_, other, _, _) => ex(s"Value $other cannot be used for an input of type InputArraySchema")
    }

    go
  }

  implicit def s: Aux[String, CommandPartBuilder] = at[String] { _ => { case (_, _, _, _) => CommandPartsList.empty } }
}
