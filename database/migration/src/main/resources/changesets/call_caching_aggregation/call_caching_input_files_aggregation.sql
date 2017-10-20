UPDATE CALL_CACHING_AGGREGATION_ENTRY ccae
  INNER JOIN (SELECT cche.CALL_CACHING_ENTRY_ID,
                UPPER(MD5(GROUP_CONCAT(
                              CONCAT(cche.HASH_KEY, cche.HASH_VALUE) ORDER BY cche.HASH_KEY ASC SEPARATOR ''
                          ))) as hash_aggregation
              FROM CALL_CACHING_HASH_ENTRY cche
              WHERE HASH_KEY LIKE 'input: File%'
              GROUP BY cche.CALL_CACHING_ENTRY_ID) cche_aggregated
    ON ccae.CALL_CACHING_ENTRY_ID = cche_aggregated.CALL_CACHING_ENTRY_ID
SET ccae.INPUT_FILES_AGGREGATION = cche_aggregated.hash_aggregation
