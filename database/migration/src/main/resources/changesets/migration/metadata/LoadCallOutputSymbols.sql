INSERT INTO TMP_SYMBOL (
  WORKFLOW_EXECUTION_UUID,
  SYMBOL_NAME,
  SYMBOL_SCOPE,
  SYMBOL_INDEX,
  SYMBOL_ATTEMPT,
  WDL_VALUE,
  WDL_TYPE
)
  SELECT
    WORKFLOW_EXECUTION_UUID,
    s.NAME,
    e.CALL_FQN,
    e.IDX,
    e.ATTEMPT,
    s.WDL_VALUE,
    s.WDL_TYPE
  FROM SYMBOL s
    JOIN WORKFLOW_EXECUTION we ON
                                 we.WORKFLOW_EXECUTION_ID = s.WORKFLOW_EXECUTION_ID
    JOIN TMP_EXECUTION_MIGRATION e ON
                            e.CALL_FQN = s.SCOPE AND
                            e.WORKFLOW_EXECUTION_ID = s.WORKFLOW_EXECUTION_ID AND
                            s.`INDEX` = e.IDX
  WHERE
    s.IO = 'OUTPUT' AND
    NOT EXISTS (                   -- filter out earlier attempts
        SELECT 1
        FROM TMP_EXECUTION_MIGRATION e3
        WHERE
          e3.WORKFLOW_EXECUTION_ID = e.WORKFLOW_EXECUTION_ID AND
          e3.CALL_FQN = e.CALL_FQN AND
          e3.IDX = e.IDX AND
          e3.ATTEMPT > e.ATTEMPT);
