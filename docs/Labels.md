Labels in Cromwell are a way to group together workflows that are related or associated to each other.
For example, if you ran workflows to analyze data from a diabetes study, you can assign the label `project:diabetes-cohort`.  

**Custom Labels JSON**

In order to assign labels to a workflow, the first step is to create a JSON file with key-value pairs that define a label. For the example above, the labels JSON should look like:

```
{
  "project":"diabetes-cohort"
}
```

When choosing key-value pairs, it's important to make sure you're adhering to Cromwell supported label syntax below.  

There are two ways to add labels to a workflow.  
Labels can be assigned to workflows upon workflow submission, by setting the `customLabels` parameter of the [submit endpoint](/api/POST_api_workflows_version) or setting the `-l` argument when running in [command line](/CommandLine) mode.

Labels can be added to existing workflows by using the [labels patch endpoint](/api/PATCH_api_workflows_version_id_labels/).

After adding labels to your workflows, you can take advantage of features like [query](/api/GET_api_workflows_version_query) to filter tagged workflows. The Google backend supports labelling cloud resources and you can learn more about that [here](/backends/Google/#google-labels).

#### Label Format

When labels are supplied to Cromwell, it will fail any request containing invalid label strings. Below are the requirements for a valid label key/value pair in Cromwell:

* Label keys and values can't contain characters other than `[a-z]`, `[0-9]` or `-`.
* Label keys must start with `[a-z]` and end with `[a-z]` or `[0-9]`.
* Label values must start and end with `[a-z]` or `[0-9]`.
* Label keys may not be empty but label values may be empty.
* Label key and values have a max char limit of 63.

Google has a different schema for label syntax requirements, where label key and value strings must match the regex `[a-z]([-a-z0-9]*[a-z0-9])?` and be no more than 63 characters in length.
For [default labels](/backends/Google/#google-labels) applied by the Google backend, Cromwell will modify workflow/task/call names to fit the schema, according to the following rules:

* Any capital letters are converted to lowercase.
* Any character which is not one of `[a-z]`, `[0-9]` or `-` will be replaced with `-`.
* If the start character does not match `[a-z]` then prefix with `x--`
* If the final character does not match `[a-z0-9]` then suffix with `--x`
* If the string is too long, only take the first 30 and last 30 characters and add `---` between them.
