SELECT
    SUM( ( (machine_type.num_cpus * preemptible.multiplier * (0.033174 / 3600) * TIMESTAMP_DIFF(end_ts.METADATA_TIMESTAMP, start_ts.METADATA_TIMESTAMP, SECOND) ) + (machine_type.memory_mib * preemptible.multiplier * (0.004446 / 3600) / 1024.0 * TIMESTAMP_DIFF(end_ts.METADATA_TIMESTAMP, start_ts.METADATA_TIMESTAMP, SECOND) ) )) AS total_cost
FROM (
         SELECT
             *,
             CAST(SPLIT(METADATA_VALUE,'-')[
                 OFFSET
                     (1)] AS INT64) AS num_cpus,
             CAST(SPLIT(METADATA_VALUE,'-')[
                 OFFSET
                     (2)] AS INT64) AS memory_mib
         FROM
             `broad-dsde-prod-analytics-dev.warehouse_dev.cromwell_metadata_2022`
         WHERE
                 METADATA_KEY = 'jes:machineType'
           AND METADATA_TIMESTAMP < '2022-01-02') machine_type
         JOIN (
    SELECT
        *
    FROM
        `broad-dsde-prod-analytics-dev.warehouse_dev.cromwell_metadata_2022`
    WHERE
            METADATA_KEY = 'start'
      AND METADATA_TIMESTAMP < '2022-01-02') start_ts
              ON
                  (machine_type.WORKFLOW_EXECUTION_UUID = start_ts.WORKFLOW_EXECUTION_UUID
                      AND machine_type.CALL_FQN = start_ts.CALL_FQN
                      AND machine_type.JOB_SCATTER_INDEX = start_ts.JOB_SCATTER_INDEX
                      AND machine_type.JOB_RETRY_ATTEMPT = start_ts.JOB_RETRY_ATTEMPT)
         JOIN (
    SELECT
        *
    FROM
        `broad-dsde-prod-analytics-dev.warehouse_dev.cromwell_metadata_2022`
    WHERE
            METADATA_KEY = 'end'
      AND METADATA_TIMESTAMP < '2022-01-02') end_ts
              ON
                  ( machine_type.WORKFLOW_EXECUTION_UUID = end_ts.WORKFLOW_EXECUTION_UUID
                      AND machine_type.CALL_FQN = end_ts.CALL_FQN
                      AND machine_type.JOB_SCATTER_INDEX = end_ts.JOB_SCATTER_INDEX
                      AND machine_type.JOB_RETRY_ATTEMPT = end_ts.JOB_RETRY_ATTEMPT)
         JOIN (
    SELECT
        *,
        (CASE METADATA_VALUE
             WHEN 'false' THEN 1
             WHEN 'true' THEN 0.1
            END
            ) AS multiplier
    FROM
        `broad-dsde-prod-analytics-dev.warehouse_dev.cromwell_metadata_2022`
    WHERE
            METADATA_KEY = 'preemptible'
      AND METADATA_TIMESTAMP < '2022-01-02') preemptible
              ON
                  ( machine_type.WORKFLOW_EXECUTION_UUID = preemptible.WORKFLOW_EXECUTION_UUID
                      AND machine_type.CALL_FQN = preemptible.CALL_FQN
                      AND machine_type.JOB_SCATTER_INDEX = preemptible.JOB_SCATTER_INDEX
                      AND machine_type.JOB_RETRY_ATTEMPT = preemptible.JOB_RETRY_ATTEMPT)
