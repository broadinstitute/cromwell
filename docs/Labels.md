Every call run on the JES backend is given certain labels by default, so that Google resources can be queried by these labels later. The current default label set automatically applied is:

| Key | Value | Example | Notes |
|-----|-------|---------|-------|
| cromwell-workflow-id | The Cromwell ID given to the root workflow (i.e. the ID returned by Cromwell on submission) | cromwell-d4b412c5-bf3d-4169-91b0-1b635ce47a26 | To fit the required [format](#label-format), we prefix with 'cromwell-' |
| cromwell-sub-workflow-name | The name of this job's sub-workflow | my-sub-workflow | Only present if the task is called in a subworkflow. |
| wdl-task-name | The name of the WDL task | my-task | |
| wdl-call-alias | The alias of the WDL call that created this job | my-task-1 | Only present if the task was called with an alias. |

## Custom Labels File

Custom labels can also be applied to every call in a workflow by specifying a custom labels file when the workflow is submitted. This file should be in JSON format and contain a set of fields: `"label-key": "label-value" `. For example:
```
{
  "label-key-1": "label-value-1",
  "label-key-2": "label-value-2",
  "label-key-3": "label-value-3"
}
```

## Label Format

When labels are supplied to Cromwell, it will fail any request containing invalid label strings. Below are the requirements for a valid label key/value pair in Cromwell:
- Label keys and values can't contain characters other than `[a-z]`, `[0-9]` or `-`.
- Label keys must start with `[a-z]` and end with `[a-z]` or `[0-9]`.
- Label values must start and end with `[a-z]` or `[0-9]`.
- Label keys may not be empty but label values may be empty.
- Label key and values have a max char limit of 63.

Google has a different schema for labels, where label key and value strings must match the regex `[a-z]([-a-z0-9]*[a-z0-9])?` and be no more than 63 characters in length.
For automatically applied labels, Cromwell will modify workflow/task/call names to fit the schema, according to the following rules:
- Any capital letters are lowercased.
- Any character which is not one of `[a-z]`, `[0-9]` or `-` will be replaced with `-`.
- If the start character does not match `[a-z]` then prefix with `x--`
- If the final character does not match `[a-z0-9]` then suffix with `--x`
- If the string is too long, only take the first 30 and last 30 characters and add `---` between them.
