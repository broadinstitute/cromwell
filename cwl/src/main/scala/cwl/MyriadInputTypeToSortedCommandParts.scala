package cwl

import cwl.CommandLineTool.{CommandBindingSortingKey, CommandPartsList, SortKeyAndCommandPart, StringOrInt}
import cwl.CwlType.CwlType
import cwl.SchemaDefRequirement.SchemaDefTypes
import cwl.command.ParentName
import shapeless.{Coproduct, Poly1}
import wom.values.{WomArray, WomCoproductValue, WomObjectLike, WomOptionalValue, WomValue}
import cats.syntax.option._
import cats.syntax.traverse._
import cats.instances.list._
import cats.instances.option._

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
  type CommandPartBuilder = (Option[InputCommandLineBinding], WomValue, CommandBindingSortingKey, Boolean, ExpressionLib, SchemaDefRequirement) => CommandPartsList

  implicit def m: Aux[MyriadInputInnerType, CommandPartBuilder] = {
    at {
      innerType => {
        (binding, womValue, sortingKey, hasShellCommandRequirement, expressionLib, schemaDefRequirement) =>
          innerType
            .fold(MyriadInputInnerTypeToSortedCommandParts)
            .apply(binding, womValue, sortingKey, hasShellCommandRequirement, expressionLib, schemaDefRequirement).fold(
              throw new RuntimeException(s"inner type to command part failed!  on type $innerType w/ womValue $womValue"))(identity)
      }
    }
  }

  implicit def am: Aux[Array[MyriadInputInnerType], CommandPartBuilder] = {
    at {
      types => {
        (binding, womValue, sortingKey, hasShellCommandRequirement, expressionLib, schemaDefRequirement) =>
          def lookupTypes(innerTypes: Array[MyriadInputInnerType]) = {
            innerTypes.toList.map{
              _.fold(MyriadInputInnerTypeToSortedCommandParts)
                .apply(binding, womValue, sortingKey, hasShellCommandRequirement, expressionLib, schemaDefRequirement)
            }.reduce(_ orElse _)
          }
          val result: Option[CommandPartsList] =
            types.partition(_.select[CwlType].contains(CwlType.Null)) match {
              case (Array(_), types) => lookupTypes(types) orElse Some(List.empty[SortKeyAndCommandPart])
              case (Array(), types) => lookupTypes(types)
            }
         result.fold(
            throw new RuntimeException(s"could not produce command line parts from input $womValue and types $types"))(
            identity
          )
      }
    }
  }
}

object MyriadInputInnerTypeToSortedCommandParts extends Poly1 {

  type CommandPartBuilder = (Option[InputCommandLineBinding], WomValue, CommandBindingSortingKey, Boolean, ExpressionLib, SchemaDefRequirement) => Option[CommandPartsList]

  import Case._

  // Primitive type: we just need to create a command part from the binding if there's one here.
  implicit def ct: Aux[CwlType, CommandPartBuilder] = {
    at {
      cwlType => {
        case (_, WomOptionalValue(_, None), _, _, _, _) => List.empty.some
        case (inputBinding, womValue, key, hasShellCommandRequirement, expressionLib, _) => {
          (cwlType, womValue) match {
            case (CwlType.Null, _) => List.empty.some
            case (_, WomOptionalValue(_, None)) => List.empty.some
            case (_,_) => inputBinding.toList.map(_.toCommandPart(key, womValue, hasShellCommandRequirement, expressionLib)).some
          }
        }
      }
    }
  }

  implicit def irs: Aux[InputRecordSchema, CommandPartBuilder] = at[InputRecordSchema] { irs => {
    def go: CommandPartBuilder = {

      //If the value is optional and is supplied, recurse over the value provided
      case (inputBinding, WomCoproductValue(_, value), sortingKey, hasShellCommandRequirement, expressionLib, schemaDefRequirement) =>
        go(inputBinding, value, sortingKey, hasShellCommandRequirement, expressionLib, schemaDefRequirement)

      //If the value is optional and is supplied, recurse over the value provided
      case (inputBinding, WomOptionalValue(_, Some(value)), sortingKey, hasShellCommandRequirement, expressionLib, schemaDefRequirement) =>
        go(inputBinding, value, sortingKey, hasShellCommandRequirement, expressionLib, schemaDefRequirement)

      // If it's optional and there's no value, do nothing
      case (_, WomOptionalValue(_, None), _, _, _, _) => List.empty.some

      // If there's no input binding and no input bindings within the irs, do nothing
      case (None, _, _, _, _, _) if !irs.fields.exists(_.exists(_.inputBinding.isDefined)) => List.empty.some

      case (inputBinding, objectLike: WomObjectLike, sortingKey, hasShellCommandRequirement, expressionLib, schemaDefRequirement) =>
        // If there's an input binding, make a SortKeyAndCommandPart for it
        val sortingKeyFromInputBindingFromInputParameter: Option[SortKeyAndCommandPart] =
          inputBinding.map(_.toCommandPart(sortingKey, objectLike, hasShellCommandRequirement, expressionLib))

        // iterate through the fields and fold over their type
        val partsFromFields:Option[CommandPartsList] =
          irs.fields.toList.flatten.
            flatTraverse{
              case InputRecordField(name, tpe, _, inputBinding, _) =>
                // Parse the name to get a clean id
                val parsedName = FullyQualifiedName(name)(ParentName.empty).id

                // The field name needs to be added to the key after the input binding (as per the spec)
                // Also start from the key from the input binding if there was one
                val fieldSortingKey =
                  sortingKeyFromInputBindingFromInputParameter.
                    map(_.sortingKey).
                    getOrElse(sortingKey).
                    append(inputBinding, Coproduct[StringOrInt](parsedName))

                val innerValueOption: Option[WomValue] = objectLike.values.get(parsedName)

                val folded = innerValueOption.map(tpe.fold(MyriadInputTypeToSortedCommandParts).
                  apply(inputBinding, _, fieldSortingKey.asNewKey, hasShellCommandRequirement, expressionLib, schemaDefRequirement))
                folded
            }

        sortingKeyFromInputBindingFromInputParameter.fold(
          partsFromFields )(
          sortingKeyFromInputBindingFromInputParameter =>
            partsFromFields.map{k => List(sortingKeyFromInputBindingFromInputParameter.copy(nestedCommandParts = k)) }
        )
      case (_, other, _, _, _, _) => throw new RuntimeException(s"Value $other cannot be used for an input of type InputRecordSchema")
    }

    go
  }}

  implicit def ies: Aux[InputEnumSchema, CommandPartBuilder] = at[InputEnumSchema] { inputEnumSchema =>
    def go: CommandPartBuilder = {
      case (_, WomOptionalValue(_, None), _, _, _, _) => List.empty.some
      case (inputBinding, value, sortingKey, hasShellCommandRequirement, expressionLib, _) =>
        // If there's an input binding, make a SortKeyAndCommandPart for it
        val fromInputBinding =
          inputBinding.map(_.toCommandPart(sortingKey, value, hasShellCommandRequirement, expressionLib)).toList

        val fromIes = inputEnumSchema.inputBinding.map(_.toCommandPart(sortingKey, value, hasShellCommandRequirement, expressionLib)).toList

        (fromInputBinding ++ fromIes).some
    }
    go
  }

  implicit def ias: Aux[InputArraySchema, CommandPartBuilder] = at[InputArraySchema] { ias =>
    def go: CommandPartBuilder = {
      case (inputBindingFromInputParameterParent, WomOptionalValue(_, Some(value)), sortingKey, hasShellCommandRequirement, expressionLib, schemaDefRequirement) =>
        go(inputBindingFromInputParameterParent, value, sortingKey, hasShellCommandRequirement, expressionLib, schemaDefRequirement)

      //If it's optional and there's no value, do nothing
      case (_, WomOptionalValue(_, None), _, _, _, _) => List.empty.some

      // If there's no input binding and no input bindings for the ias, do nothing
      case (None, _, _, _, _, _) if ias.inputBinding.isEmpty => List.empty.some

      case (inputBinding, WomArray.WomArrayLike(womArray: WomArray), sortingKey, hasShellCommandRequirement, expressionLib, schemaDefRequirement) =>


        // If there's an input binding, make a SortKeyAndCommandPart for it
        val sortKeyFromInputBindingFromInputerParameterParent: Option[SortKeyAndCommandPart] =
          inputBinding.
            map(_.toCommandPart(sortingKey, womArray, hasShellCommandRequirement, expressionLib))

        // Now depending on whether we have an itemSeparator and/or valueFrom or not, we're going to recurse over each element of the array (or not).
        // See http://www.commonwl.org/v1.0/CommandLineTool.html#CommandLineBinding
        if (inputBinding.flatMap(_.itemSeparator).isDefined || inputBinding.flatMap(_.valueFrom).isDefined) {
          // If there's an item separator or a valueFrom we can stop here.
          // When the command part is instantiated (see CommandLineBindingCommandPart) it will evaluate the valueFrom (if defined) and join the items together (if there's an itemSeparator).
          sortKeyFromInputBindingFromInputerParameterParent.toList.some
        } else {
          // If neither valueFrom nor itemSeparator were defined, we need to process each item of the array
          val fromArray: Option[CommandPartsList] = womArray.value.zipWithIndex.toList.flatTraverse({
            case (item, index) =>
              // Update the sorting key with the binding position (if any), add the index
              val itemSortingKey = {

              // The index needs to be added to the key after the input binding (as per the spec)
              // Also start from the key from the input binding if there was one
                val sortKey = sortKeyFromInputBindingFromInputerParameterParent.fold(sortingKey)(_.sortingKey)

                sortKey.append(ias.inputBinding, Coproduct[StringOrInt](index))
              }

              // Even if the item doesn't have an explicit input binding, it should appear in the command so create a default empty one
              //there is an explicit input binding!
              val arrayItemInputBinding =
                ias.
                  inputBinding.
                  orElse(Option(InputCommandLineBinding.default))

              // Fold over the item type of each array element
              Option(ias.items.fold(MyriadInputTypeToSortedCommandParts).apply(arrayItemInputBinding, item, itemSortingKey.asNewKey, hasShellCommandRequirement, expressionLib, schemaDefRequirement))
          })
          sortKeyFromInputBindingFromInputerParameterParent.fold(fromArray){
            sortKeyFromInputBindingFromInputerParameterParent =>
              fromArray.map{k => List(sortKeyFromInputBindingFromInputerParameterParent.copy(nestedCommandParts = k)) }
          }
        }

      case (_, other, _, _, _, _) => throw new RuntimeException(s"Value $other cannot be used for an input of type InputArraySchema")
    }

    go
  }

  implicit def s: Aux[String, CommandPartBuilder] = {
    at {
      s => {
        case (_, WomOptionalValue(_, None), _, _, _, _) => List.empty.some
        case (commandLineBindingFromInput, womValue, sortingKey, boolean, expressionLib, schemaDefRequirement) => {
          val womType: SchemaDefTypes = schemaDefRequirement.lookupCwlType(s).getOrElse(throw new RuntimeException(s"Looked for type $s in custom types $schemaDefRequirement but no match was found!"))

          womType.fold(this).apply(commandLineBindingFromInput, womValue, sortingKey, boolean, expressionLib, schemaDefRequirement)
        }
      }
    }
  }
}
