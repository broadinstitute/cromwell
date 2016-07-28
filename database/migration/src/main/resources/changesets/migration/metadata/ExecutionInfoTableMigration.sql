INSERT INTO METADATA_JOURNAL (
  WORKFLOW_EXECUTION_UUID,
  METADATA_KEY,
  CALL_FQN,
  JOB_SCATTER_INDEX,
  JOB_RETRY_ATTEMPT,
  METADATA_VALUE,
  METADATA_VALUE_TYPE,
  METADATA_TIMESTAMP
)
  SELECT
    WORKFLOW_EXECUTION_UUID,
    CASE ei.INFO_KEY
    WHEN '$log_stdout' THEN 'stdout'
    WHEN '$log_stderr' THEN 'stderr'
    WHEN 'JES_RUN_ID' THEN 'jobId'
    WHEN 'JES_STATUS' THEN 'backendStatus'
    WHEN 'SGE_JOB_NUMBER' THEN 'jobNumber'
    WHEN 'LSF_JOB_NUMBER' THEN 'jobNumber'
    ELSE
      IF(ei.INFO_KEY LIKE '$log_%',                       -- backend log
         CONCAT(
             'backendLogs:',                                     -- prepend metadata prefix
             SUBSTRING(ei.INFO_KEY, 6, LENGTH(ei.INFO_KEY) - 5)  -- remove $log_ prefix from the key
         ),
         ei.INFO_KEY -- Just put the key otherwise
      ) END,
    CALL_FQN,
    IDX,
    ATTEMPT,
    ei.INFO_VALUE,
    'string',
    NOW()
  FROM EXECUTION_INFO ei
    LEFT JOIN TMP_EXECUTION_MIGRATION e ON ei.EXECUTION_ID = e.EXECUTION_ID
    JOIN WORKFLOW_EXECUTION we ON we.WORKFLOW_EXECUTION_ID = e.WORKFLOW_EXECUTION_ID
  WHERE
    ei.INFO_VALUE IS NOT NULL;
