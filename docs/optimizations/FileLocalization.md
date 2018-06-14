# The 'localization_optional' Optimization

## Scope

The 'localization_optional' optimization applies to `File` and `File?` inputs on `task`s. It allows
you to save time and money by identifying `File` inputs which do not need to be localized for the task to succeed.

## Condition

The optimization signals to Cromwell that a task has been written in such a way that:

 * The task **will work** if Cromwell does localize the specified file input 
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
  }
  
  parameter_meta {
    foo_file: {
      localization_optional: true
    }
  }
  
  command <<<
    # This tool must work for **BOTH** local file paths **AND** object store URL values:
    java -jar my_tool.jar ~{foo_file}
  >>>
}
```

## Backend Support

This optimization is currently only applied to localization in the Pipelines API (GCE) backends.
