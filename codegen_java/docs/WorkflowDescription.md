

# WorkflowDescription


## Properties

| Name | Type | Description | Notes |
|------------ | ------------- | ------------- | -------------|
|**valid** | **Boolean** | Indicates that the workflow is valid and that the inputs, if provided, are compatible with the workflow. |  |
|**errors** | **List&lt;String&gt;** | The set of validation failure messages |  |
|**validWorkflow** | **Boolean** | Indicates whether the workflow file is valid by itself. If inputs are provided, they are not considered when calculating this field; if inputs are not provided, the value is identical to &#x60;valid&#x60;. |  |
|**name** | **String** | For a source file with one workflow and zero or more tasks, the name of the workflow. For a single task, the name of the task. For a source file with multiple tasks but no workflows, the empty string. |  |
|**inputs** | [**List&lt;ToolInputParameter&gt;**](ToolInputParameter.md) | A list of inputs for this tool |  |
|**outputs** | [**List&lt;ToolOutputParameter&gt;**](ToolOutputParameter.md) | A list of outputs for this tool |  |
|**submittedDescriptorType** | [**DescriptorTypeAndVersion**](DescriptorTypeAndVersion.md) |  |  |
|**isRunnableWorkflow** | **Boolean** | Indicates whether this file can be run on its own (e.g. a WDL workflow) |  |
|**importedDescriptorTypes** | **List&lt;String&gt;** |  |  [optional] |
|**meta** | **Object** |  |  [optional] |
|**parameterMeta** | **Object** |  |  [optional] |



