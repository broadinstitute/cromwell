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
    s.SCOPE,
    e.IDX,
    e.ATTEMPT,
    s.WDL_VALUE,
    s.WDL_TYPE
  FROM SYMBOL s
    JOIN WORKFLOW_EXECUTION we ON
                                 we.WORKFLOW_EXECUTION_ID = s.WORKFLOW_EXECUTION_ID
    LEFT JOIN TMP_EXECUTION_MIGRATION e ON
                            e.CALL_FQN = s.SCOPE AND
                            e.WORKFLOW_EXECUTION_ID = s.WORKFLOW_EXECUTION_ID
                            -- Don't join on index here because inputs have index null (-1), but we want to duplicate them as many times as there are indices for a given execution
  WHERE
    s.IO = 'INPUT';
