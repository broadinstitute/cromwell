

# ValueType

The type expected for a given value.

## Properties

| Name | Type | Description | Notes |
|------------ | ------------- | ------------- | -------------|
|**typeName** | [**TypeNameEnum**](#TypeNameEnum) | The type of this value |  [optional] |
|**optionalType** | [**ValueType**](ValueType.md) |  |  [optional] |
|**arrayType** | [**ValueType**](ValueType.md) |  |  [optional] |
|**mapType** | [**MapValueType**](MapValueType.md) |  |  [optional] |
|**tupleTypes** | [**List&lt;ValueType&gt;**](ValueType.md) |  |  [optional] |
|**objectFieldTypes** | [**List&lt;ValueTypeObjectFieldTypesInner&gt;**](ValueTypeObjectFieldTypesInner.md) |  |  [optional] |



## Enum: TypeNameEnum

| Name | Value |
|---- | -----|
| STRING | &quot;String&quot; |
| FILE | &quot;File&quot; |
| DIRECTORY | &quot;Directory&quot; |
| FLOAT | &quot;Float&quot; |
| INT | &quot;Int&quot; |
| BOOLEAN | &quot;Boolean&quot; |
| OPTIONAL | &quot;Optional&quot; |
| ARRAY | &quot;Array&quot; |
| TUPLE | &quot;Tuple&quot; |
| MAP | &quot;Map&quot; |
| OBJECT | &quot;Object&quot; |
| PAIR | &quot;Pair&quot; |



