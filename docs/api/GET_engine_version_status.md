Cromwell will track the status of underlying systems on which it depends, typically connectivity to the database and to Docker Hub. The `status` endpoint will return the current status of these systems. 

Response:

This endpoint will return an `Internal Server Error` if any systems are marked as failing or `OK` otherwise. The exact response will vary based on what is being monitored but adheres to the following pattern. Each system has a boolean `ok` field and an optional array field named `messages` which will be populated with any known errors if the `ok` status is `false`:

```
{
  "System Name 1": {
    "ok": false,
    "messages": [
      "Unknown status"
    ]
  },
  "System Name 2": {
    "ok": true
  }
}

```