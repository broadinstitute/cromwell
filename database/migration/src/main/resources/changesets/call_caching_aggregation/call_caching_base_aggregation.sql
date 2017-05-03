INSERT INTO CALL_CACHING_AGGREGATION_ENTRY(
  CALL_CACHING_ENTRY_ID,
  BASE_AGGREGATION
)
  SELECT cche.CALL_CACHING_ENTRY_ID,
  UPPER(MD5(GROUP_CONCAT(
                cche.HASH_VALUE ORDER BY cche.HASH_KEY ASC SEPARATOR ''
            ))) as aggregated
FROM CALL_CACHING_HASH_ENTRY cche
WHERE HASH_KEY NOT LIKE 'input: File%'
AND
      (HASH_KEY NOT LIKE 'runtime attribute: %' OR
       HASH_KEY = 'runtime attribute: continueOnReturnCode' OR
       HASH_KEY = 'runtime attribute: docker' OR
       HASH_KEY = 'runtime attribute: failOnStderr'
      )
        
GROUP BY cche.CALL_CACHING_ENTRY_ID;
