version 1.1

# Prints more lines than the default Cloud Logging cap (10,000) to both stdout and stderr, then checks that:
#  - the stdout/stderr files still contain every line, since the cap only limits what is forwarded to Cloud Logging;
#  - Cloud Logging received each stream's 10,000th line and its truncation note, but not its 10,001st line.
workflow cloud_logging_line_cap {
  call spam

  call check_cloud_logging { input: marker = spam.marker, spam_stdout = spam.stdout_file, spam_stderr = spam.stderr_file }

  output {
    Int stdout_lines = check_cloud_logging.stdout_lines
    Int stderr_lines = check_cloud_logging.stderr_lines
    Boolean cloud_logging_stdout_line_10000 = check_cloud_logging.stdout_line_10000
    Boolean cloud_logging_stderr_line_10000 = check_cloud_logging.stderr_line_10000
    Boolean cloud_logging_stdout_note = check_cloud_logging.stdout_note
    Boolean cloud_logging_stderr_note = check_cloud_logging.stderr_note
    Boolean cloud_logging_stdout_line_10001 = check_cloud_logging.stdout_line_10001
    Boolean cloud_logging_stderr_line_10001 = check_cloud_logging.stderr_line_10001
  }
}

task spam {
  # Every line carries a per-run marker so the Cloud Logging check finds only this run's lines.
  command <<<
    cat /proc/sys/kernel/random/uuid > marker.txt
    marker=$(cat marker.txt)
    seq 1 20000 | sed "s/^/out-$marker-/"
    seq 1 20000 | sed "s/^/err-$marker-/" >&2
  >>>

  output {
    String marker = read_string("marker.txt")
    # Too big for read_lines (128,000-byte limit), so check_cloud_logging counts the lines instead.
    File stdout_file = stdout()
    File stderr_file = stderr()
  }

  runtime {
    docker: "mirror.gcr.io/almalinux:latest"
    predefinedMachineType: "c3d-standard-4"
    disks: "local-disk 50 SSD"
  }

  meta {
    # A cache hit would produce no new log entries for the check to find.
    volatile: true
  }
}

task check_cloud_logging {
  input {
    String marker
    File spam_stdout
    File spam_stderr
  }

  # Never print the marker here: this task's own output also goes to Cloud Logging and would match the queries.
  command <<<
    set -euo pipefail
    shopt -s inherit_errexit
    marker="~{marker}"

    wc -l < "~{spam_stdout}" > stdout_lines.txt
    wc -l < "~{spam_stderr}" > stderr_lines.txt

    # Prints the Batch job UID of the spam task's line number $2 on stream $1, or nothing if it isn't in Cloud Logging.
    line_job() {
      gcloud logging read "logName:\"batch_task_logs\" AND textPayload:\"$1-$marker-$2\"" \
        --freshness=1d --limit=1 --format='value(labels.job_uid)'
    }

    # Prints true if the truncation note for stream $1 reached Cloud Logging from Batch job $2.
    # Added Logs Viewer permission to `centaur@broad-dsde-cromwell-dev` enable
    has_note() {
      local found
      found=$(gcloud logging read "logName:\"batch_task_logs\" AND labels.job_uid=\"$2\" AND textPayload:\"[Cromwell] $1 truncated after\"" \
        --freshness=1d --limit=1 --format='value(insertId)')
      if [ -n "$found" ]; then echo true; else echo false; fi
    }

    # The spam task finished before this task started, so its logs are normally already there. Retry briefly in
    # case ingestion is lagging; this only adds time when a check would otherwise fail.
    for attempt in 1 2 3 4 5; do
      stdout_job=$(line_job out 10000)
      stderr_job=$(line_job err 10000)
      stdout_note=false
      stderr_note=false
      if [ -n "$stdout_job" ]; then stdout_note=$(has_note stdout "$stdout_job"); fi
      if [ -n "$stderr_job" ]; then stderr_note=$(has_note stderr "$stderr_job"); fi
      echo "Attempt $attempt: stdout_note=$stdout_note stderr_note=$stderr_note"
      if [ "$stdout_note" = true ] && [ "$stderr_note" = true ]; then break; fi
      sleep 30
    done

    if [ -n "$stdout_job" ]; then echo true; else echo false; fi > stdout_line_10000.txt
    if [ -n "$stderr_job" ]; then echo true; else echo false; fi > stderr_line_10000.txt
    echo "$stdout_note" > stdout_note.txt
    echo "$stderr_note" > stderr_note.txt

    # Each stream's note is emitted in place of its line 10,001, so once the note has arrived, an uncapped
    # line 10,001 would have too.
    stdout_extra_job=$(line_job out 10001)
    stderr_extra_job=$(line_job err 10001)
    if [ -n "$stdout_extra_job" ]; then echo true; else echo false; fi > stdout_line_10001.txt
    if [ -n "$stderr_extra_job" ]; then echo true; else echo false; fi > stderr_line_10001.txt
  >>>

  output {
    Int stdout_lines = read_int("stdout_lines.txt")
    Int stderr_lines = read_int("stderr_lines.txt")
    Boolean stdout_line_10000 = read_boolean("stdout_line_10000.txt")
    Boolean stderr_line_10000 = read_boolean("stderr_line_10000.txt")
    Boolean stdout_note = read_boolean("stdout_note.txt")
    Boolean stderr_note = read_boolean("stderr_note.txt")
    Boolean stdout_line_10001 = read_boolean("stdout_line_10001.txt")
    Boolean stderr_line_10001 = read_boolean("stderr_line_10001.txt")
  }

  runtime {
    docker: "mirror.gcr.io/google/cloud-sdk:stable"
    predefinedMachineType: "c3d-standard-4"
    disks: "local-disk 50 SSD"
  }

  meta {
    volatile: true
  }
}
