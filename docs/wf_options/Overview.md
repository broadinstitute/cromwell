_For the Doc-A-Thon_
**Questions to answer and things to consider:**

1. Who is visiting the workflow options page?
*They may have been directed by another page in the docs*
2. What do they need to know first?

3. Is all the important information there? If not, add it!

4. Are there things that don't need to be there? Remove them.

5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---


When running a workflow from the [command line](/commandline#run) or [REST API](/restapi#post-apiworkflowsversion), one may specify a JSON file that toggles various options for running the workflow.  From the command line, the workflow options is passed in as the third positional parameter to the 'run' subcommand.  From the REST API, it's an optional part in the multi-part POST request.  See the respective sections for more details.

Example workflow options file:

```json
{
  "jes_gcs_root": "gs://my-bucket/workflows",
  "google_project": "my_google_project",
  "refresh_token": "1/Fjf8gfJr5fdfNf9dk26fdn23FDm4x"
}
```


