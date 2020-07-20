# Workflow jobs: retry decision logic

* **Word count:** 132

## Key concepts

* Workflow may consist of multiple jobs.
* Job "retry" actually means creating a new job.
* Job may have attributes which are being stored in the `JOB_KEY_VALUE_ENTRY` table in the database, identified by a 
`ScopedKey`, which comprises workflow id, call fully qualified name, job index, job attempt number, and attribute name.
* There are 2 types of retries:
  * backend-specific retries (e.g., VM preemption in PAPI)
  * general retries 
* Backend-specific and general retries have separate retry counters, which are being stored in `JOB_KEY_VALUE_ENTRY` in
the end of the job execution attempt, and pre-fetched from the table in the beginning of the next attempt.

The retry logic is shown on the sequence diagram below with the example of PAPIv2 backend and VM preemption as an example 
of backend-specific retry reason.

![Job retry logic (example with PAPIv2 and VM preemption)](Workflow_job_retry_logic_(example_with_PAPIv2_and_VM_preemption).png)
