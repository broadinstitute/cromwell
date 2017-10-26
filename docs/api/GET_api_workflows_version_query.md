_For the Doc-A-Thon_  
**Questions to answer and things to consider:**

1. Who is visiting the API/Query page?  

2. What do they need to know first?  

3. Is all the important information there? If not, add it!  
*Why would someone want to know about queries? What can it help them do?*
4. Are there things that don't need to be there? Remove them.  

5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---


This endpoint allows for querying workflows based on the following criteria:

* `name`
* `id`
* `status`
* `label`
* `start` (start datetime with mandatory offset)
* `end` (end datetime with mandatory offset)
* `page` (page of results)
* `pagesize` (# of results per page)

Names, ids, and statuses can be given multiple times to include
workflows with any of the specified names, ids, or statuses. When
multiple names are specified, any workflow matching one of the names
will be returned. The same is true for multiple ids or statuses. When
more than one label is specified, only workflows associated to all of
the given labels will be returned. 

When a combination of criteria are specified, for example querying by 
names and statuses, the results must return workflows that match one of 
the specified names and one of the statuses. Using page and pagesize will
enable server side pagination.

Valid statuses are `Submitted`, `Running`, `Aborting`, `Aborted`, `Failed`, and `Succeeded`.  `start` and `end` should
be in [ISO8601 datetime](https://en.wikipedia.org/wiki/ISO_8601) format with *mandatory offset* and `start` cannot be after `end`.

cURL:

```
$ curl "http://localhost:8000/api/workflows/v1/query?start=2015-11-01T00%3A00%3A00-04%3A00&end=2015-11-04T00%3A00%3A00-04%3A00&status=Failed&status=Succeeded&page=1&pagesize=10"
```

HTTPie:

```
$ http "http://localhost:8000/api/workflows/v1/query?start=2015-11-01T00%3A00%3A00-04%3A00&end=2015-11-04T00%3A00%3A00-04%3A00&status=Failed&status=Succeeded&page=1&pagesize=10"
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 133
Content-Type: application/json; charset=UTF-8
Date: Tue, 02 Jun 2015 18:06:56 GMT
Server: spray-can/1.3.3

{
  "results": [
    {
      "name": "w",
      "id": "fdfa8482-e870-4528-b639-73514b0469b2",
      "status": "Succeeded",
      "end": "2015-11-01T07:45:52.000-05:00",
      "start": "2015-11-01T07:38:57.000-05:00"
    },
    {
      "name": "hello",
      "id": "e69895b1-42ed-40e1-b42d-888532c49a0f",
      "status": "Succeeded",
      "end": "2015-11-01T07:45:30.000-05:00",
      "start": "2015-11-01T07:38:58.000-05:00"
    },
    {
      "name": "crasher",
      "id": "ed44cce4-d21b-4c42-b76d-9d145e4d3607",
      "status": "Failed",
      "end": "2015-11-01T07:45:44.000-05:00",
      "start": "2015-11-01T07:38:59.000-05:00"
    }
  ]
}
```

Labels have to be queried in key and value pairs separated by a colon, i.e. `label-key:label-value`. For example, if a batch of workflows was submitted with the following labels JSON:
```
{
  "label-key-1" : "label-value-1",
  "label-key-2" : "label-value-2"
}
```

A request to query for succeeded workflows with both labels would be:

cURL:
```
$ curl "http://localhost:8000/api/workflows/v1/query?status=Succeeded&label=label-key-1:label-value-1&label=label-key-2:label-value-2
```

HTTPie:
```
$ http "http://localhost:8000/api/workflows/v1/query?status=Succeeded&label=label-key-1:label-value-1&label=label-key-2:label-value-2
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 608
Content-Type: application/json; charset=UTF-8
Date: Tue, 9 May 2017 20:24:33 GMT
Server: spray-can/1.3.3

{
    "results": [
        {
            "end": "2017-05-09T16:07:30.515-04:00", 
            "id": "83fc23d5-48d1-456e-997a-087e55cd2e06", 
            "name": "wf_hello", 
            "start": "2017-05-09T16:01:51.940-04:00", 
            "status": "Succeeded"
        }, 
        {
            "end": "2017-05-09T16:07:13.174-04:00", 
            "id": "7620a5c6-a5c6-466c-994b-dd8dca917b9b", 
            "name": "wf_goodbye", 
            "start": "2017-05-09T16:01:51.939-04:00", 
            "status": "Succeeded"
        }
    ]
}
```

Query data is refreshed from raw data periodically according to the configuration value `services.MetadataService.metadata-summary-refresh-interval`.
This interval represents the duration between the end of one summary refresh sweep and the beginning of the next sweep.  If not specified the
refresh interval will default to 2 seconds.  To turn off metadata summary refresh, specify an infinite refresh interval value with "Inf".

```
services {
  MetadataService {
    config {
      metadata-summary-refresh-interval = "10 seconds"
    }
  }
}
```