# The 'localization_optional' Optimization

Available in Cromwell version 33 and higher.

## Scope

The 'localization_optional' optimization can be applied to a task's individual input declarations containing files, specifically `File` and `File?` values and any complex types containing them. 
It allows you to save time and money by identifying files which do not need to be localized for the task to succeed.

## Condition

The optimization signals to Cromwell that a task has been written in such a way that:
 
 * The task **will work** if Cromwell does localize the specified file inputs
   * For example if a file is localized for a local dockerized execution environment.

**And**:

 * The task will **also** work if Cromwell **does not** localize the same file input
   * For example the file remains in a cloud object store and the command is constructed using its URL rather than a local path.

## Effect on File Localization

If the [backend](#backend-support) has been set up to respect `localization_optional`, Cromwell will 
choose not to localize the appropriate file input.

### Effect on Call Caching:

None! 

Files marked for optional localization are still treated in exactly the same way as other `File` inputs for call caching.

## Language Support

### WDL 1.0 (or later)

In a WDL 1.0 `task`, this optimization is specified by adding a `localization_optional` field to 
an input's entry in the task's `parameter_meta` section. Here's an example:

```wdl
task nio_task {
  input {
    File foo_file
    File bar_file
  }
  
  parameter_meta {
    foo_file: {
      description: "a foo file",
      localization_optional: true
    }
    bar_file: {
      description: "a bar file"
    }
  }
  
  command <<<
    # This tool must work for **BOTH** local file paths **AND** object store URL values:
    java -jar my_tool_1.jar ~{foo_file}
    
    # Because the optimization is not applied to 'bar_file' in parameter_meta, this file **WILL** be localized:
    java -jar my_tool_2.jar ~{bar_file}
  >>>
}
```

## Backend Support

This optimization is currently only applied to localization in the Pipelines API (GCE) backends.
