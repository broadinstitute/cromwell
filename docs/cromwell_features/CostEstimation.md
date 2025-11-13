Cromwell provides estimated cloud costs for running and terminal workflows. This estimate describes costs that have already been incurred by this workflow and does not attempt to predict future costs.

## Usage

This feature is currently available for **GCP only**. To enable it, add this line to your configuration file:
```hocon
services.GcpCostCatalogService.config.enabled = true
```

Ensure that you have [enabled the Cloud Billing API](https://cloud.google.com/billing/v1/how-tos/catalog-api) for your project and that your Google credentials are available in your environment via the `GOOGLE_APPLICATION_CREDENTIALS` env var or [another method](https://cloud.google.com/docs/authentication/provide-credentials-adc).

Workflows run with this functionality enabled will produce non-zero cost when the `/api/workflows/{version}/{workflowId}/cost` endpoint is called. For example:
```
/api/workflows/v1/00001111-2222-3333-aaaa-bbbbccccdddd/cost

{
  "cost": 0.045,
  "currency": "USD",
  "errors": [
    "string"
  ],
  "id": "00001111-2222-3333-aaaa-bbbbccccdddd",
  "status": "Running"
}
```
 * `cost` (number) the current estimation of cloud costs incurred by this workflow (see caveats below). This includes cores and memory for this workflow and all its subworkflows.
 * `currency` (string) the currency the cost is displayed in, currently always USD
 * `errors` (list of strings) errors encountered when computing cost, if any. For example, inability to look up cloud costs for a task via the GCP API.
 * `id` (string) the UUID of the workflow
 * `status` (string) the current status of the workflow

## Caveats

The initial implementation of this feature focuses on the most common cost centers and use case: VM cores and memory priced according to [the public pricing catalog](https://cloud.google.com/billing/v1/how-tos/catalog-api). The following are NOT currently reflected in Cromwell's estimate:
 * Any cloud discounts that you or your organization may have negotiated with GCP
 * GPU costs
 * Disk costs
 * Data egress costs

## Implementation

When each task starts running in GCP, Cromwell notes:
 * What region the compute is in
 * Whether the machine is a spot instance or dedicated compute
 * What type the machine is (N1, N2, N2D, etc)
 * How many cores the machine has
 * How much memory the machine has

This information is used to determine the cost per hour of the VM. That hourly cost, along with the time the VM started running and the time the VM stopped running, are saved in the Cromwell metadata table. When the cost API endpoint is called, Cromwell gathers this cost, start, stop metadata for each task in the workflow, recursing down to all subworkflows. (For tasks that are in progress, the current wallclock time is used as the end time.) This metadata is used to compute the current estimated cost for each task; those are summed to return the estimated cost for the root workflow.
