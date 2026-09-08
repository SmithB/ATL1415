# MAAP platform documentation -- notes for the ATL14/15 transition

Read 2026-09-08 against <https://docs.maap-project.org/en/latest/> (the "latest" build;
the site banner says MAAP is transitioning from the ADE to the MAAP Hub, and a separate
`last-ade-release` build exists for the legacy ADE).  Latest release note on the site is
5.1.0, 2026-03-03.

Conventions used below:
  STATE   -- the docs say this, with the URL.
  SOURCE  -- read out of the installed `maap-py` 4.2.0 at
             /srv/conda/envs/notebook/lib/python3.13/site-packages/maap/ , not from the docs.
             Included where the docs are silent and the library is authoritative.
  INFER   -- my deduction, not stated anywhere.
  SILENT  -- the docs do not cover this.  Recorded so nobody re-searches.
  CONFLICT -- the docs disagree with something we established by experiment, or with
             themselves.

---

## 0.  Findings that change something for this project

Most consequential first.

### 0.1  The public worker queues have NO time limit; the private org queue we asked for has 24 h
STATE.  <https://docs.maap-project.org/en/latest/system_reference_guide/dps_queues.html>
gives, verbatim:

| Queue | Memory | Time limit | Guest access | AWS Instance Type |
|---|---|---|---|---|
| maap-dps-sandbox | 8gb | 10 min | yes | t3.large or t3a.large |
| maap-dps-worker-8gb | 8gb | unlimited | no | t3.large or t3a.large |
| maap-dps-worker-16gb | 16gb | unlimited | no | t3.xlarge or t3a.xlarge |
| maap-dps-worker-32gb | 32gb | unlimited | no | r5.xlarge or r5a.xlarge |

and for private organization queues: "24 hr time limit on individual jobs (can be increased
upon request)", no throttling, and "custom resources can be requested".

Why it matters: the transition doc has been treating walltime as "a property of the queue"
(Q7) without knowing the numbers.  The numbers are: public = unlimited, org = 24 h.  A 4-hour
SLURM tile fits either.  The sandbox's 10 min is why smoke tests must be trivial.

### 0.2  The instance types pin the core counts -- `maap-dps-worker-32gb` is 4 vCPU, not more
STATE (instance types, same URL).  INFER (the vCPU/RAM per instance type, from AWS's published
specs, not from MAAP docs): t3.large / t3a.large = 2 vCPU / 8 GiB; t3.xlarge / t3a.xlarge =
4 vCPU / 16 GiB; r5.xlarge / r5a.xlarge = 4 vCPU / 32 GiB.

This corroborates the measured "sandbox worker reports 2 cores" (t3.large) and predicts that
`run.sh`'s `threads=$(nproc)` will be **4** on our registered default queue
`maap-dps-worker-32gb` -- the same as the 4 tasks the discover SLURM runs use.  No change
needed to run.sh; this just removes the guess from algorithm_config.yml's "Sizing is still a
guess" comment.

Note the t3/t3a queues are *burstable* instances (INFER): sustained CPU on a t3.xlarge is
throttled to a baseline once credits are exhausted.  A long, CPU-bound solve on the 16gb queue
may therefore run much slower than 4 full cores.  The 32gb queue (r5) is not burstable.

### 0.3  DPS workers are spot instances; a fan-out of thousands MUST have a requeue path
STATE. <https://docs.maap-project.org/en/latest/troubleshooting/dps_error_codes.html> --
the entire error-code table is one row:

| error code | explanation | resolution |
|---|---|---|
| exit code 143 | The system lost access to a resource due to spot-instance interruption | Re-run the job(s) later or with a different resource |

INFER: at a fan-out of thousands of multi-hour tiles, exit 143 will be a routine, non-zero
fraction of the failures and is *not* a code bug.  `check_MAAP_jobs.py` (Q10) should special-case
exit code 143 as "requeue, do not triage", separately from real exceptions.  DPS does not
retry automatically as far as the docs say (SILENT on retries).

### 0.4  The platform ceiling is 4,000 concurrent jobs, and per-member job limits are an org setting
STATE. <https://docs.maap-project.org/en/latest/getting_started/running_at_scale.html>:
"MAAP is configured to run up to 4,000 concurrent jobs."
STATE. <https://docs.maap-project.org/en/latest/system_reference_guide/organizations.html>:
Org Maintainers "can add/remove members and set per-member job limits."

So the Q11 "max jobs in flight" policy has two real constraints: a 4,000 platform-wide ceiling
(shared with every other MAAP user) and a per-member limit that an Org Maintainer of the
`icesat-2` org sets.  Ask the org maintainer what our per-member limit is -- that is the number
`--max_in_flight` has to respect, and it is not visible from `getQueues()`.

### 0.5  There is NO batch/bulk submission API.  One HTTP POST per job, full stop.
SOURCE (`maap/dps/DpsHelper.py:submit_job`, `maap/dps/execute.xml`): `submitJob` renders a
single OGC WPS 2.0 `<wps:Execute ... mode="sync">` document with
`<ows:Identifier>job-{algo_id}:{version}</ows:Identifier>` and one `<wps:Input id="NAME">`
carrying a `<wps:LiteralValue><![CDATA[value]]></wps:LiteralValue>` per kwarg, then POSTs it.
Every kwarg that is not `algo_id`/`version`/`inputs` becomes an algorithm input; `username` is
injected if absent (defaulting to `anonymous`).
STATE. The docs' only "batch" advice is
<https://docs.maap-project.org/en/latest/getting_started/running_at_scale.html>: "You can press
Submit Job repeatedly", and the DPS tutorial's "Use maap.py library to loop through datasets
(e.g., 1000 jobs from 1000 input files)".

SILENT: job arrays, bulk submit, submit-from-file, async submission.  The Q11 submitter must be
a Python loop with its own rate limiting; nothing in the platform does it for us.
INFER: because the Execute request is `mode="sync"`, each submitJob blocks on an HTTP round
trip.  A thousands-job fan-out is latency-bound; thread the submitter or accept ~1 job/second.

### 0.6  `getJobMetrics` returns peak memory and duration -- use it to size the queue empirically
SOURCE (`maap/dps/dps_job.py`, `DPSJob` properties populated by `retrieve_metrics()` /
`maap.getJobMetrics(jobid)`):
`machine_type`, `architecture`, `machine_memory_size`, `directory_size`, `operating_system`,
`job_start_time`, `job_end_time`, `job_duration_seconds`, `cpu_usage`, `cache_usage`,
`mem_usage`, `max_mem_usage`, `swap_usage`, `read_io_stats`, `write_io_stats`, `sync_io_stats`,
`async_io_stats`, `total_io_stats`, `outputs`, `response_code`, `error_details`.
STATE: the docs mention `getJobMetrics` only in passing; the only memory guidance they give is
"profile your own script first" --
<https://docs.maap-project.org/en/latest/technical_tutorials/user_data/memory-profiling-python.html>:
"This is useful when you have working code and you want to estimate the size of the DPS worker
to be used."

`max_mem_usage` and `job_duration_seconds` are exactly what Q4/Q18 and the queue-sizing TBD
need, measured on the real thing.  Record both in the ledger after each job completes.

### 0.7  The documented output prefix contains the ALGORITHM VERSION, not the job tag
CONFLICT (internal to the docs, and with our Q9 assumption).
- STATE. <https://docs.maap-project.org/en/latest/system_reference_guide/jobs_maappy.html>:
  "The output data will be put into a folder named for your `algo_id` and the `identifier`."
- STATE. DPS tutorial FAQ,
  <https://docs.maap-project.org/en/latest/technical_tutorials/dps_tutorial/dps_tutorial_demo.html>:
  outputs are under `/home/jovyan/my-private-bucket/dps_output` organised "by algorithm name,
  job tag, and a set of folders organized by date and time".
- STATE. But the *worked example path* in the Jobs UI page,
  <https://docs.maap-project.org/en/latest/system_reference_guide/jobsui.html>, is
  `/projects/my-private-bucket/dps_output/run-dps-test_ubuntu/delay10/2023/05/10/15/13/27/250000`
  and in that page's own submission example `algo_id="run-dps-test_ubuntu"`,
  `version="delay10"`, `identifier="test-job"`.  So the second component is the **version**,
  and `test-job` appears nowhere in the path.

INFER: `dps_output/<algo_id>/<algorithm_version>/YYYY/MM/DD/HH/MM/SS/<microseconds>/`, i.e.
the same for every tile of a given registration, with only a microsecond-resolution timestamp
distinguishing jobs.  Either way it is not addressable from the tile coordinates, which is the
premise of Q9.  **The Q9 decision (have the job write its own deterministic
`s3://maap-ops-workspace/ben_smith/ATL14_processing/rel<NNN>/<period>/<hemi>/<REGION>/{prelim,matched}/E%d_N%d.h5`
key) stands and is confirmed necessary.**  Do not try to reconstruct DPS's own prefix.
The authoritative per-job answer is `maap.getJobResult(job_id)`, which returns the product URLs
(HTTP, S3 and AWS-console forms) -- STATE, jobs_maappy.html.

### 0.8  Nothing in the docs supports a job localizing another job's outputs
SILENT.  There is no mechanism for "give me the outputs of job X as an input to job Y", no
job dependency/DAG support, no workflow engine.  File inputs are URLs resolved at submission
time only, and "if the value provided [for] a parameter marked as file during registration is
not a valid url, DPS will report an error" (DPS tutorial FAQ).
This confirms the Q8 answer (a): the `matched` job must fetch its 9 prelim tiles itself with
`aws s3 cp` from the deterministic Q9 prefix.  There is no platform feature we are missing.

### 0.9  `/tmp` exists on DPS workers and is the fastest disk; home/`/projects` does not
STATE. <https://docs.maap-project.org/en/latest/system_reference_guide/disk_guide.html>:
- Home directory: `/projects/` on the ADE (EFS), `/home/jovyan/` on the Hub (EBS on NFS),
  "Backed up Daily, backups for 30 days", moderate performance, **not available in DPS jobs**,
  "should be used for code, small sample files and other documents".
- Local disk `/tmp`: "Fastest performance", "There is no backup", "cleared on every reboot of
  your workspace", and "Also exists on DPS workers so you can make it part of your scripts
  reliably".  "The disk space for /tmp is not infinite" -- no number given.
- Buckets: "Available in both Workspaces and DPS", "Slower initial read (higher latency)",
  "High throughput (many parallel reads)", "we keep old versions for 30 days (this includes
  deletions)", and **"Use of S3 paths, direct reads are encouraged for best performance"**
  (i.e. the docs themselves prefer `s3://` / `/vsis3/` over the FUSE mount).

INFER: the Q8(a) `aws s3 cp` of the 9 prelim tiles should land in `/tmp` or the job dir, not
anywhere under a mount, and `disk_space: 20GB` in algorithm_config.yml is what reserves it.

### 0.10  `GDAL_DISABLE_READDIR_ON_OPEN=EMPTY_DIR` is the documented S3 read speedup
STATE. <https://docs.maap-project.org/en/latest/technical_tutorials/access/accessing_cod.html>:
"speed up GDAL reads from S3 buckets by skipping sidecar files" by setting
`GDAL_DISABLE_READDIR_ON_OPEN` to `'EMPTY_DIR'`.
Directly applicable: every `/vsis3/` mask and ancillary-grid read in a tile job pays a LIST per
open otherwise.  Candidate for an `export` in `run.sh`.  (Caveat, INFER: `EMPTY_DIR` also
suppresses discovery of genuine sidecars -- `.aux.xml`, `.msk`, and for the RGI `.db` masks
nothing is at risk, but a `.tif` with an external overview `.ovr` would lose it.)

### 0.11  Secrets exist and a running job can read them
STATE. <https://docs.maap-project.org/en/latest/system_reference_guide/jobs_maappy.html>,
"Passing Credentials for Other Services into Jobs":
```python
maap.secrets.add_secret("<SECRET_NAME>", "<SECRET_VALUE>")
maap.secrets.get_secrets()
maap.secrets.get_secret("<SECRET_NAME>")
maap.secrets.delete_secret("<SECRET_NAME>")
```
and "Inside your Algorithm code, you will use the maap-py `secrets.get_secret("SECRET_NAME")`
method" to retrieve values during job execution.  Secrets are "encrypted during transmission".
Introduced in release 4.1.0 (2024-10-02) for e.g. Google Earth Engine credentials.
INFER: we do not need this -- `earthdata_s3_credentials` already works on the worker without a
stored secret -- but it is the documented escape hatch if a future input needs a real
credential, and it would be the right place for an Earthdata token if the s3credentials route
ever stops working.
SILENT: secret scope (per-user vs per-org), expiry, size limits, whether they appear as env vars.

---

## 1.  DPS job submission at scale

STATE -- `submitJob` example, <https://docs.maap-project.org/en/latest/system_reference_guide/jobs_maappy.html>:
```python
maap.submitJob(identifier="test-job",
               algo_id="run-dps-test_ubuntu",
               version="delay10",
               queue="maap-dps-worker-8gb",
               input_file="https://raw.githubusercontent.com/MAAP-Project/dps-unit-test/main/README.md")
```
returns `{'status': 'success', 'http_status_code': 200, 'job_id': '<uuid>'}`.

SOURCE -- real signature is
`submitJob(self, identifier, algo_id, version, queue, retrieve_attributes=False, **kwargs)`.
`retrieve_attributes=True` makes an extra round trip to populate the returned `DPSJob`; leave it
False in a fan-out loop.  `kwargs['username']` is overwritten from `profile.account_info()`.
(Both consistent with what we established.)

STATE -- throttling, dps_queues.html: "The platform team will start throttling users ability to
submit a batch of jobs on public queues to about **10 jobs per hour**."  Private org queues:
"no throttling".
STATE -- concurrency ceiling: 4,000 platform-wide (running_at_scale.html).
STATE -- per-member job limits are set by Org Maintainers (organizations.html).

Requesting an organizational queue: STATE, organizations.html -- "Creating organizations,
designating maintainers, and assigning job queues to an organization remain MAAP administrator
functions", and "Only job queues assigned to your organization(s) will be allowed for use when
registering an algorithm or submitting jobs."  dps_queues.html says to "Contact the platform
team to request private queues", with the specs quoted in 0.1.
SILENT: the docs give **no email address, ticket URL or form** for that request.  The route is
whatever channel we already used (the request is already in flight per the transition doc).
Worth asking for explicitly, given "custom resources can be requested": a queue sized to the
near-pole AA tiles (more RAM and more vCPU than r5.xlarge) and a walltime above 24 h.

SILENT: job arrays; bulk submit; submission rate-limit HTTP status codes or Retry-After headers;
what the throttle does when exceeded (reject, queue, or drop); max in-flight per user as a number.

CONFLICT with our established queue list: `getQueues()` returns `maap-dps-worker-64gb` and
`maap-dps-worker-32vcpu-64gb`, which the dps_queues.html table does **not** list.  The docs
table is stale.  INFER: those two exist but their instance types, time limits and org
eligibility are undocumented -- if we get access, ask for the instance type explicitly rather
than assuming.  `-32vcpu-64gb` would be the obvious home for near-pole AA tiles.

---

## 2.  Job monitoring at scale

STATE -- `listJobs` query parameters, jobs_maappy.html:

| Parameter | Type | Values |
|---|---|---|
| algo_id | string | valid string |
| end_time | string, `2024-01-01T00:00:00.000Z` | jobs completed from that time to now |
| get_job_details | bool | False or True |
| offset | number | integer |
| page_size | number | integer |
| priority | number | 0-9 |
| queue | string | valid string |
| start_time | string | jobs started from that time to now |
| status | string | `job-queued`, `job-started`, `job-completed`, `job-failed`, `job-revoked`, `job-offline` |
| tag | string | valid string |
| version | string | valid string |

SOURCE -- the installed signature is keyword-only and defaults `get_job_details=True`,
`offset=0`, **`page_size=10`**; it raises `ValueError` unless `algo_id` and `version` are
either both given or both omitted, and internally collapses them to
`job_type=f"{algo_id}:{version}"`.  `priority` is documented on the page but is **not** in the
4.2.0 `listJobs` signature -- CONFLICT (docs list a parameter the library will reject as an
unexpected keyword).

STATE -- status mapping, jobs_maappy.html: Accepted <- job-queued, Running <- job-started,
Success <- job-completed, Failed <- job-offline/job-failed, job-revoked <- job-revoked.

Finding jobs by identifier: the `identifier=` passed to `submitJob` is the **job tag**, and
`listJobs(tag=...)` filters on it (STATE for the parameter; INFER for the identity of
identifier and tag, supported by running_at_scale.html: "Job Tags serve as identifiers for
tracking" with the example tag `test_run_2024b`).
INFER, actionable: tag every tile job `ATL1415_<region>_<step>_E<x>_N<y>` so the ledger can be
rebuilt from `listJobs` if the local CSV is lost -- cheap insurance on the Q10 design, which
otherwise depends entirely on a local file.
SILENT: whether `tag` matching is exact or a prefix/substring match.  Test it.

STATE -- other monitoring calls, jobs_maappy.html:
`maap.getJobStatus(job_id)` (returns XML), `maap.getJobResult(job_id)` (XML with output URLs in
HTTP, S3 and AWS-console forms), `maap.cancelJob(job_id)` -- "may take several minutes before
the cancellation takes effect".
SOURCE -- also `maap.getJob(jobid)` (returns a populated `DPSJob`), `job.wait_for_completion()`,
and `DPSJob.retrieve_metrics()` (see 0.6).  `dps_job.py` wraps the status/metrics polls in
`backoff` retry decorators, so transient API errors are already handled inside maap-py.

STATE -- Jobs UI, jobsui.html: sortable by queued/start/end time; "use the search bar to filter
the job list down to jobs containing the user-provided string in any of the fields shown";
per-job Outputs tab with a "Products" path and `_stderr.txt`, `_stdout.txt`, and context/dataset/met
JSON; cancel is available only for queued or running jobs (queued -> deleted, running -> stopped
and marked `job-revoked`); "Copy Jupyter Notebook Code" generates the maap-py submission call.
STATE -- release 4.1.1 (2024-10-23): "Increased Jupyter Jobs UI record limit to 200."  INFER:
the UI is unusable as the primary monitor for a thousands-job fan-out; the ledger + `listJobs`
is the right design.
STATE -- release 3.1.4 (2024-01-22) added a "dps-job-management shared workspace for job
tracking"; release 4.2.0 (2025-03-03) added "job duration tracking" to the Jobs UI.

SILENT: bulk cancel via API (only the single-job `cancelJob`); any requeue/retry API; job
priority as a submit-time argument (the `priority` field appears only as a listJobs filter);
notification/webhook on completion.

Triaged jobs: SILENT.  Nothing in the docs describes the
`s3://maap-ops-workspace/dataset/triaged_job/v1.4.0/<job dir>/` tree we found, beyond
disk_guide.html naming a `triaged-jobs` bucket "for debugging failed DPS jobs".  Our
behaviour-derived layout is the only documentation of it.

---

## 3.  Outputs

STATE -- DPS tutorial FAQ, dps_tutorial_demo.html, verbatim:
> "Since the jobs on DPS are run on a machine on the cloud, your local workspace directories are
> not available. It is important to pass any files required as inputs for your algorithm using
> the `File` parameter type. Any outputs that need to be saved should be placed in a directory
> called `output`. When a parameter is registered as a file input, DPS downloads the
> corresponding value provided by the user as a file and places it in a directory called
> `input`. It is important to note that if the value provided a parameter marked as file during
> registration is not a valid url, DPS will report an error. Note: Both `input` and `output`
> directories are relative to your run script."

> "File management i.e. files required for input and files stored as outputs on S3 are taken
> care of by the DPS. To locate the files created as an output from your job, look into the
> `/home/jovyan/my-private-bucket/dps_output` dir on your workspace and navigate to the
> algorithm type and time of run. You can also construct the output path of your files by
> looking at the job info on the Jobs UI or by running `maap.getJobResult('job_id')`"

Output prefix: see 0.7.  Observed example path (jobsui.html):
`/projects/my-private-bucket/dps_output/run-dps-test_ubuntu/delay10/2023/05/10/15/13/27/250000`.
INFER: `my-private-bucket` is the FUSE view of `s3://maap-ops-workspace/<username>/`, so the S3
form is `s3://maap-ops-workspace/<username>/dps_output/<algo_id>/<version>/<YYYY/MM/DD/HH/MM/SS/us>/`.
(Checked 2026-09-08: `s3://maap-ops-workspace/ben_smith/dps_output/` does not exist yet -- no
DPS job of ours has produced output.)

Naming/finding outputs afterwards: SILENT beyond `getJobResult`.  There is no
"name my output" field, no output-prefix override, no manifest.
Dataset/product ingestion: STATE -- the only ingestion documentation is STAC metadata guidance,
<https://docs.maap-project.org/en/latest/technical_tutorials/user_data/stac_metadata.html>,
which specifies required STAC collection-level fields (id, version, title, description,
providers, keywords, spatial/temporal extent, platforms, instruments, processing level, license)
and item-level fields (id, collection, geometry, bbox, datetime, links, assets with href and
type).  It does **not** describe the submission workflow, an API, or where the data must live.
Catalogs, <https://docs.maap-project.org/en/latest/technical_tutorials/searching.html>: MAAP
STAC <https://stac.maap-project.org>, NASA CMR STAC <https://cmr.earthdata.nasa.gov/stac/ALL>,
ESA MAAP STAC <https://catalog.maap.eo.esa.int/catalogue/>.  `cmr.maap-project.org` was
deprecated 2023-05-01.  There is also a curated bucket `nasa-maap-data-store`
(<https://docs.maap-project.org/en/latest/technical_tutorials/access/aws_access.html>).
INFER: ATL14/15 product delivery goes to NSIDC as it always has; MAAP STAC ingestion is not on
our path and should not be spent time on.

Later job localizing earlier jobs' outputs: see 0.8.  SILENT.

---

## 4.  Inputs and localization

STATE -- two registration input kinds, dps_tutorial_demo.html:
- **File inputs**: "Parameters downloaded as files and placed in `input/` directory. Must be
  valid URLs or DPS reports an error."
- **Positional inputs**: "Command-line parameters ... passed directly to run script."
- "The order of inputs in registration must match the order passed to the run script."

STATE -- run-script pattern the docs bless:
```bash
basedir=$(dirname "$(readlink -f "$0")")
mkdir -p output
conda run --live-stream --name dps_tutorial python ${basedir}/gdal_wrapper.py \
  --input_file ${INPUT_FILENAME} \
  --output_file output/${OUTPUT_FILENAME} \
  --outsize ${REDUCTION_SIZE}
```
FAQ: `basedir` exists because "it is not possible to know the absolute path of your script
before execution".  (Our `run.sh` already does the same via `readlink -f`.)
Note the docs use `conda run --live-stream`; we use `conda run --no-capture-output`.  INFER:
equivalent for our purposes -- both stream to the job's `_stdout.txt`.

CONFLICT / SILENT, and this is the big one for us: the docs describe file localization as
"DPS **downloads** the corresponding value ... as a file and places it in a directory called
`input`".  They say nothing about the shared read-only cache, nothing about **symlinks**, and
nothing about `/data/work/cache/<md5>/`.  Our established fact -- that an
`s3://maap-ops-workspace/...` file input is localized as a *symlink* into a shared cache, which
is why `find -L` is load-bearing in run.sh -- is entirely behaviour-derived and is contradicted
in spirit by the docs' "downloads".  Keep the comment block in run.sh; it is the only record.

SILENT, all of it:
- whether `s3://` vs `https://` matters for a file input (the docs' only example is an
  `https://raw.githubusercontent.com/...` URL).  Our experience says `s3://` works.
- any per-input size limit.
- any per-job limit on the number of file inputs.
- whether a directory or prefix can be given instead of a single object (it cannot, as far as
  anything documented goes -- reinforcing Q8(a)).
- whether the cache is shared between jobs, deduplicated, or evicted.
- `config` inputs.  Our algorithm_config.yml declares `config: []`; the docs never mention a
  `config` input class at all.

INFER, practical consequence: because a `file` input is one URL and there is no directory
input, the *only* documented way to give a job many files is many declared inputs, and every
one of them must be declared at registration time -- i.e. a fixed arity.  The 3x3 prelim
neighbourhood has fixed arity 9, so a 9-file-input registration is technically possible, but
it would require nine URLs computed per submission and would break at region edges where
neighbours do not exist.  `aws s3 cp` inside run.sh (Q8(a)) remains the right call.

---

## 5.  Credentials and data access

SOURCE -- `maap/AWS.py`, the complete surface of `maap.aws`:
```python
maap.aws.requester_pays_credentials(expiration=60*60*12)      # default 12 h
maap.aws.s3_signed_url(bucket, key, expiration=60*60*12)
maap.aws.earthdata_s3_credentials(endpoint_uri)               # double-url-quotes endpoint_uri
maap.aws.workspace_bucket_credentials()
```
`earthdata_s3_credentials` adds a `"DAAC"` key (the endpoint hostname) to whatever the DAAC
returns.  All four are thin GETs against MAAP member-API endpoints using the maap-py API header
(i.e. `MAAP_PGT`), which is why they work on a DPS worker.

STATE -- `workspace_bucket_credentials()`,
<https://docs.maap-project.org/en/latest/system_reference_guide/accessing_bucket_data.html>:
returns a `credentials` object with `aws_access_key_id`, `aws_secret_access_key`,
`aws_session_token`, `expires_at` (example `"2025-03-03T18:00:00Z"`), plus an
**authorized S3 paths array** whose entries have `bucket`, `prefix`, `uri`, `type`
(`workspace` or `org`) and `access` (`read_write` or `read_only`).  Example URIs in the docs:
`s3://maap-ops-workspace/maap_user`, `s3://shared-project-bucket/team-data`,
`s3://public-reference-data/smap/v9`.  Usage:
```python
import boto3
creds = resp["credentials"]
session = boto3.Session(
    aws_access_key_id=creds["aws_access_key_id"],
    aws_secret_access_key=creds["aws_secret_access_key"],
    aws_session_token=creds["aws_session_token"],
)
s3 = session.client("s3")
```
Release 4.1.1 (2024-10-23) "expanded S3 permissions for `aws.workspace_bucket_credentials()`".
INFER: this is the supported way for a DPS job to get *write* credentials to
`s3://maap-ops-workspace/ben_smith/...` if the worker's own instance role turns out not to
allow the Q9 deterministic-prefix write.  Worth knowing before the first production fan-out --
the smoke test should include one `aws s3 cp` write to the target prefix.
SILENT: the docs "[do] not explicitly specify how DPS jobs should utilize these credentials",
and give no credential lifetime other than the example timestamp.

STATE -- reading DAAC data.  There are **three different documented routes**, and none of them
is the one we settled on:
1. `maap.searchGranule(...)` + `maap.getData()` / `maap.downloadGranule(online_access_url, ...)`
   over HTTPS, requiring the user to have authorised the DAAC application in their Earthdata
   Login profile --
   <https://docs.maap-project.org/en/latest/technical_tutorials/access/accessing_external_data.html>,
   <https://docs.maap-project.org/en/latest/technical_tutorials/access/accessing_data.html>.
2. The `maap-data-reader` assumed role, via an SSM parameter --
   <https://docs.maap-project.org/en/latest/technical_tutorials/access/direct_access.html>:
```python
def assume_role_credentials(ssm_parameter_name):
    session = boto3.Session()
    ssm = session.client('ssm', "us-west-2")
    parameter = ssm.get_parameter(Name=ssm_parameter_name, WithDecryption=True)
    sts = session.client('sts')
    assumed_role_object = sts.assume_role(
        RoleArn=parameter['Parameter']['Value'], RoleSessionName='TutorialSession')
    return assumed_role_object['Credentials']
```
   with `ssm_parameter_name = "/iam/maap-data-reader"`, and then
   `fsspec.filesystem("s3", key=..., secret=..., token=..., requester_pays=requester_pays)` or
   `rasterio.session.AWSSession(..., requester_pays=requester_pays)`.  Covers GES DISC, LPDAAC,
   NSIDC, ORNL, PO.DAAC.  Added in release 4.0.0 (2024-07-03).
3. Ambient Hub credentials with `os.environ["AWS_REQUEST_PAYER"] = "requester"` --
   <https://docs.maap-project.org/en/latest/technical_tutorials/access/external_access_from_hub.html>,
   explicitly "an experimental feature".

CONFLICT with our established fact.  Routes 2 and 3 both carry a hard caveat -- "This tutorial
must be run within MAAP's ADE" and "This tutorial must be run within MAAP's hub to assume the
necessary permissions" -- and route 2 says **GES DISC, LPDAAC and NSIDC are Requester Pays
buckets** requiring `requester_pays=True`.  We established that
`maap.aws.earthdata_s3_credentials('https://data.nsidc.earthdatacloud.nasa.gov/s3credentials')`
works from a DPS worker and returns `accessKeyId`/`secretAccessKey`/`sessionToken`.
INFER: these are different mechanisms and both are true.  Route 2 assumes a MAAP-owned IAM role
and reads the DAAC's *public/requester-pays* S3 bucket, which is why it is region- and
platform-bound and bills the requester.  `earthdata_s3_credentials` proxies the DAAC's own
`/s3credentials` endpoint using the MAAP API token, returning DAAC-issued short-lived keys for
the protected bucket, and needs only `MAAP_PGT` -- hence it works on a worker.  **Our route is
the better one and the docs simply do not describe it**; the tutorials' "must be run within the
ADE/Hub" does not apply to us.  Nothing here should change our implementation.
SILENT: `earthdata_s3_credentials` is not documented on any page I read -- it exists only in the
library.  Its credential lifetime is whatever NSIDC issues (typically 1 h; INFER, not stated) --
a multi-hour tile job must be able to **refresh** them mid-run.  Flag this as a real risk for
the long AA tiles: if the solve reads ATL11 granules more than ~an hour after job start with
credentials captured at start, it will 403.

STATE -- `MAAP_PGT`,
<https://docs.maap-project.org/en/latest/system_reference_guide/personal_access_tokens.html>:
"Outside the Hub, `maap-py` reads your token from the `MAAP_PGT` environment variable":
```python
os.environ["MAAP_PGT"] = "<your-personal-access-token>"
```
PATs are created in the MAAP Console under Profile > Personal Access Tokens, with a Platform
(NASA or ESA) and an Expiration of "`1 day`, `7 days`, `30 days`, `90 days`, `1 year`, or
`No expiration`".  "Create a new token whenever one expires -- expired tokens cannot be
renewed."  SILENT on DPS workers, but this is exactly the variable we found set in the
container, and it explains why maap-py works there.
SOURCE -- `MAAP(maap_host=os.getenv('MAAP_API_HOST', 'api.maap-project.org'))`, which is why
`MAAP_API_HOST` is also set in the container.

STATE -- region: "We do need to set the default AWS region to `us-west-2`"
(aws_access.html), and the SSM/STS calls in direct_access.html are pinned to `us-west-2`.
INFER: MAAP DPS runs in us-west-2, the same region as the NSIDC Earthdata Cloud bucket and
`s3://pytmd`, so our reads are in-region and free of egress.  Not stated anywhere.

---

## 6.  Algorithm registration

STATE -- registration entry points:
`maap.register_algorithm_from_yaml_file("/home/jovyan/<algorithm_config>.yml").text`
(dps_tutorial_demo.html), or the Register Algorithm UI, which "automatically generates YAML
files in the home directory after initial registration, which can be reused for updates"
(<https://docs.maap-project.org/en/latest/system_reference_guide/algorithm_registration.html>;
the generated files land in an `algorithm-configs` folder).
`maap.deleteAlgorithm("<algorithm_name>:<branch>")` unregisters --
<https://docs.maap-project.org/en/latest/system_reference_guide/faq/delete_algorithm_from_mas.html>.
SOURCE -- also `maap.listAlgorithms()`, `maap.describeAlgorithm(algoid)`,
`maap.publishAlgorithm(algoid)` (undocumented), `maap.registerAlgorithm(dict_or_json)`.

**The docs do not publish an algorithm_config.yml schema.**  This is the single biggest
documentation gap for us.  What exists:

STATE -- the field list sketched in dps_tutorial_demo.html:
```yaml
algorithm_name:
algorithm_description:
repository_url:
repository_branch:
build_command:
run_command:
disk_space:
resource:
container_url:
file_inputs:
  - name: ...
    description: ...
    default: ...
    required: ...
positional_inputs:
  - name: ...
    description: ...
    default: ...
    required: ...
```
with "Input order in YAML must match the order of arguments passed to the run script."

CONFLICT: that field list does not match ours.  We use `algorithm_version` (not
`repository_branch`), `docker_container_url` (not `container_url`), `queue` (not `resource`),
and `inputs: {positional:, file:, config:}` (not top-level `file_inputs`/`positional_inputs`).
The MAAP-Project/dps-unit-test reference repo
(<https://raw.githubusercontent.com/MAAP-Project/dps-unit-test/main/algorithm_config.yaml>) is
older still -- `algo_name`, `version`, `environment`, `description`, `docker_url`, `inputs`
with only `name` and `download` -- and is annotated "THIS CONFIG IS AUTO-GENERATED BY ADE UI".

SOURCE, and this explains the mess: `maap/utils/algorithm_utils.py` is
```python
def read_yaml_file(algo_yaml):
    with open(algo_yaml) as fr: algo_config = yaml_load(fr, Loader=Loader)
    return validate_algorithm_config(algo_config)

def validate_algorithm_config(algo_config):
    return algo_config
```
-- i.e. `register_algorithm_from_yaml_file` performs **no validation whatsoever** and POSTs the
YAML verbatim as JSON to the MAS register endpoint.  The separate
`register_algorithm_from_yaml_file_backwards_compatible` translates the *old* keys
(`algo_name`->`algorithm_name`, `version`->`code_version`, `environment`->`environment_name`,
`description`->`algorithm_description`, `docker_url`->`docker_container_url`,
`run_command`->`script_command`, `repository_url`->`repo_url`, `inputs`->`algorithm_params`
flattened to `{field, download}`), which pins the *new* names as the canonical ones -- and
`docker_container_url` is among them, confirming our file uses the current schema and the
tutorial's `container_url` is the stale one.
INFER: the authoritative schema is the MAS server's, not the client's.  An unknown field is
silently forwarded, so a typo in algorithm_config.yml fails at registration or, worse, is
ignored.  `describeAlgorithm(algoid)` after registering is the only way to confirm what was
actually stored -- which is what our register_algorithm.py already does.

SILENT, every one of these -- **the docs contain no algorithm_config.yml field for**:
- CPU or memory requests (only `queue`/`resource`, i.e. pick a queue).
- job timeout or walltime (walltime is a queue property; 0.1).
- retries / max attempts.
- environment variables passed into the job.
- secrets bound to the algorithm (secrets are fetched at runtime by name; 0.11).
- build caching or build timeout.
- an output-prefix or output-naming field.
So there is nothing we are leaving unused.  `disk_space` (as `"20GB"`, `"5GB"`, `"100GB"`,
`"20MB"`, `"10KB"`) and `queue` are the only resource knobs.

Container images: STATE, release 4.1.0 (2024-10-02) "Introduced `maap_base` minimal container
for DPS with faster registration", and 4.2.0 (2025-03-03) "Added `maap_base` container option
as default for algorithm registration (fastest option)".  The Register Algorithm UI offers "a
dropdown selection from `maap_base` or workspace container".
SILENT: <https://docs.maap-project.org/en/latest/system_reference_guide/custom-environments.html>
covers Hub workspace environments only and says nothing about DPS base images, their registry
URLs/tags, `DOCKERIMAGE_PATH_DEFAULT`, or whether a custom-built image can be used for DPS
registration.  It points at <https://github.com/MAAP-Project/maap-workspaces> and
<https://docs.openveda.cloud/user-guide/scientific-computing/custom-environments.html>.
So our `mas.maap-project.org/root/maap-workspaces/custom_images/maap_base:v6.0.0` choice, and
the reasoning for it, has no documentation to check against.
STATE -- release 4.2.0 (2025-03-03): "Switched all images to mini-forge; conda installs now
pull exclusively from conda-forge", consistent with maap_base being Miniforge-based and with
build-env.sh's conda-forge assumption.

Build logs: SILENT.  Nothing describes where a build log lives or how to fetch it
programmatically -- consistent with our finding that build output is browser-only.

---

## 7.  Storage

STATE -- disk_guide.html, quoted in 0.9.  Adding:
- getting_started.html "MAAP Storage Options": `~/` mounted at `/home/jovyan` is "local (to
  Jupyter) file system; generally faster and more reliable", "Use this for code-related items,
  smaller data storage"; `~/my-private-bucket` is an S3 bucket, "persistent storage, but
  accessible only to you and others in a shared workspace", "for large data storage";
  `~/my-public-bucket` is "equivalent to `~/shared-buckets/<my_username>/`" and read-only to
  other users at `~/shared-buckets/<their_username>`.
- ADE->Hub migration,
  <https://docs.maap-project.org/en/latest/system_reference_guide/faq/ade_to_hub.html>: "your
  **buckets** will remain the same"; files in `my-private-bucket` / `my-public-bucket` are
  "visible from the ADE and Hub"; **default home directory quota is 150 GB** (adjustable case
  by case); default idle workspace timeout is one hour; instance size is chosen at each
  workspace launch.
- share_data.html: presigned URLs via `Command Palette -> User -> Get Presigned S3 Url` or
  right-click "Get Presigned S3 Url"; "The link will expire after 12 hours".

CONFLICT on the workspace-bucket prefix.  Three inconsistent statements:
- aws_access.html says private is `s3://maap-ops-workspace/private/<username>/...` and shared is
  `s3://maap-ops-workspace/shared/<username>/...`.
- accessing_bucket_data.html's example `authorized_s3_paths` entry is
  `s3://maap-ops-workspace/maap_user` (no `private/`).
- Observed 2026-09-08 by `aws s3 ls s3://maap-ops-workspace/`: the bucket's top level is a flat
  list of usernames (`ben_smith/`, `aimeeb/`, ...), i.e. `s3://maap-ops-workspace/<username>/`
  with **no** `private/` segment.  Our established `s3://maap-ops-workspace/ben_smith/...` is
  correct and aws_access.html is wrong or describes something else.
INFER: `~/my-private-bucket` == `s3://maap-ops-workspace/<username>/`.  A `shared/` prefix may
also exist for the public-bucket view; not verified.

mountpoint-s3 vs direct S3: STATE -- disk_guide.html, "Use of S3 paths, direct reads are
encouraged for best performance", with buckets characterised as "Slower initial read (higher
latency)" / "High throughput (many parallel reads)".  This endorses what MEMORY.md already
records (no rename on the FUSE mount; reorganize with `aws s3 mv`) and what MAAP_dps.txt does
(every static input as an `s3://` URI).  SILENT on mountpoint-s3 specifically, and on its
limitations (no rename, no partial write) -- that remains behaviour-derived.

Costs and egress: **SILENT**.  Nothing anywhere on the site about data-transfer cost, egress,
who pays for compute, per-user billing, or budget limits, other than: the requester-pays flag
for GES DISC/LPDAAC/NSIDC direct-bucket access (direct_access.html), and the implication of a
per-member job limit set by Org Maintainers.  If cost is a planning input for a thousands-job
fan-out, it has to come from the platform team, not the docs.

Versioning/backups: STATE -- buckets keep "old versions for 30 days (this includes deletions)";
home is "Backed up Daily, backups for 30 days"; `/tmp` has "no backup".

---

## 8.  Efficiency, cost, container practice

STATE:
- `GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR'` to skip sidecar probing on S3 reads
  (accessing_cod.html).  See 0.10.
- COG overviews: "Overviews are versions of the data with lower resolution, and can thus
  increase performance in applications" (accessing_cod.html).
- `maap_base` is "the fastest option" for algorithm registration (release notes 4.2.0).
- `/tmp` is "Fastest performance" and exists on DPS workers (disk_guide.html); buckets have
  "High throughput (many parallel reads)".
- Memory-profile before choosing a worker size (memory-profiling-python.html).

SILENT: build caching between registrations (INFER, and worth stating plainly: since DPS clones
`repository_url` at `algorithm_version` and builds a container per registration, every
re-registration re-runs `build-env.sh` end to end -- the ~conda+SuiteSparse cost is paid again.
Nothing documents a cache, so minimise re-registrations and pin `algorithm_version` to a tag,
not a moving branch, once production starts).
SILENT: any guidance on avoiding redundant downloads across jobs, on the shared input cache, on
in-region reads, or on request-pattern tuning for S3.

---

## 9.  Consolidated list of things the docs do NOT cover

Recorded so nobody searches for them again.

1. Batch/bulk job submission, job arrays, job dependencies or any DAG/workflow feature.
2. Any way for one job to consume another job's outputs.
3. A deterministic or overridable output prefix; output naming.
4. The input localization *mechanism* -- shared cache, symlinks, `/data/work/cache/<md5>/`.
   The docs say "downloads ... as a file", which is not what we observed.
5. Per-job limits on input count or input size; whether `s3://` is a supported file-input scheme.
6. `config` inputs (our algorithm_config.yml has `config: []`).
7. The complete algorithm_config.yml schema; the current field names anywhere in one place.
8. Any algorithm_config field for cpu, memory, timeout, retries, env vars, or secrets.
9. `maap.aws.earthdata_s3_credentials` -- not documented on any page; library-only.
10. Credential lifetimes for `earthdata_s3_credentials` and `workspace_bucket_credentials`,
    and how to refresh them inside a long-running job.
11. The `maap-dps-worker-64gb` and `maap-dps-worker-32vcpu-64gb` queues that `getQueues()`
    reports; vCPU counts per queue.
12. A named contact, form or ticket URL for requesting an organizational queue.
13. Max in-flight jobs per user; what the ~10 jobs/hr throttle does when exceeded.
14. DPS error codes other than 143; the `triaged_job` bucket layout; retry semantics.
15. Bulk cancel; requeue; job priority at submit time.
16. Where DPS build logs live or any API to fetch them.
17. DPS base-image URLs/tags, `DOCKERIMAGE_PATH_DEFAULT`, custom images for DPS registration.
18. Build caching; container layer reuse across registrations.
19. Costs, billing, egress, or data-transfer charges of any kind.
20. The AWS region MAAP DPS runs in (only inferable from `us-west-2` in the code examples).
21. Secret scope, expiry, or size limits.
22. Whether `listJobs(tag=...)` matches exactly or by substring.
23. The `/tmp` size on a DPS worker.

---

## 10.  Where the docs contradict our established facts

| # | Established (behaviour-derived) | Docs say | Assessment |
|---|---|---|---|
| 1 | A `file` input is localized as a **symlink** into `/data/work/cache/<md5>/` | "DPS downloads the corresponding value ... as a file and places it in a directory called `input`" (dps_tutorial_demo.html FAQ) | Docs are simplified, not wrong about the end state.  Our `find -L` fix stands.  Nothing to change. |
| 2 | Output prefix contains the job **tag** | jobsui.html's worked path contains the **version** (`run-dps-test_ubuntu/delay10/...`), while the same site's prose says "algo_id and the identifier" and "algorithm name, job tag" | Docs contradict themselves.  Either way the prefix is not addressable from tile coordinates -- Q9's self-written deterministic key is required. |
| 3 | `earthdata_s3_credentials` works from a DPS worker | The two DAAC-access tutorials both say "must be run within MAAP's ADE" / "within MAAP's hub" | Different mechanism (SSM assume-role vs DAAC `/s3credentials` proxy).  Our route is undocumented but real.  No change. |
| 4 | NSIDC ATL11 read works without requester-pays | direct_access.html: GES DISC, LPDAAC and NSIDC are Requester Pays, need `requester_pays=True` | Applies to the `maap-data-reader` role reading the DAAC's requester-pays bucket, not to DAAC-issued credentials for the protected bucket.  No change, but if a 403 ever appears with "requester pays" in it, this is the explanation. |
| 5 | `getQueues()` lists `maap-dps-worker-64gb` and `-32vcpu-64gb` | dps_queues.html table lists only sandbox/8gb/16gb/32gb | Docs table is stale.  Do not assume the two extra queues have documented properties. |
| 6 | Our schema: `algorithm_version`, `docker_container_url`, `queue`, `inputs:{positional,file,config}` | Tutorial sketch: `repository_branch`, `container_url`, `resource`, `file_inputs`/`positional_inputs`; dps-unit-test repo: `algo_name`, `docker_url`, `version` | Three generations of schema in the docs.  maap-py's backwards-compat key map confirms ours is current.  Ignore the tutorial's field names. |
| 7 | Workspace prefix is `s3://maap-ops-workspace/ben_smith/...` | aws_access.html: `s3://maap-ops-workspace/private/<username>/` | Verified against the live bucket 2026-09-08: top level is flat usernames.  aws_access.html is wrong. |
| 8 | `listJobs` filters | jobs_maappy.html documents a `priority` (0-9) filter | Not present in maap-py 4.2.0's `listJobs` signature; passing it raises TypeError. |

---

## 11.  Source index

Pages read in full for these notes:

- <https://docs.maap-project.org/en/latest/>
- <https://docs.maap-project.org/en/latest/getting_started.html>
- <https://docs.maap-project.org/en/latest/getting_started/getting_started.html>
- <https://docs.maap-project.org/en/latest/getting_started/maap_overview.html>
- <https://docs.maap-project.org/en/latest/getting_started/writing_code.html>
- <https://docs.maap-project.org/en/latest/getting_started/running_at_scale.html>
- <https://docs.maap-project.org/en/latest/technical_tutorials.html>
- <https://docs.maap-project.org/en/latest/technical_tutorials/dps_tutorial/dps_tutorial_demo.html>
- <https://docs.maap-project.org/en/latest/technical_tutorials/searching.html>
- <https://docs.maap-project.org/en/latest/technical_tutorials/access/accessing_data.html>
- <https://docs.maap-project.org/en/latest/technical_tutorials/access/accessing_external_data.html>
- <https://docs.maap-project.org/en/latest/technical_tutorials/access/accessing_cod.html>
- <https://docs.maap-project.org/en/latest/technical_tutorials/access/aws_access.html>
- <https://docs.maap-project.org/en/latest/technical_tutorials/access/direct_access.html>
- <https://docs.maap-project.org/en/latest/technical_tutorials/access/external_access_from_hub.html>
- <https://docs.maap-project.org/en/latest/technical_tutorials/user_data/stac_metadata.html>
- <https://docs.maap-project.org/en/latest/technical_tutorials/user_data/memory-profiling-python.html>
- <https://docs.maap-project.org/en/latest/system_reference.html>
- <https://docs.maap-project.org/en/latest/system_reference_guide/disk_guide.html>
- <https://docs.maap-project.org/en/latest/system_reference_guide/share_data.html>
- <https://docs.maap-project.org/en/latest/system_reference_guide/algorithm_registration.html>
- <https://docs.maap-project.org/en/latest/system_reference_guide/jobsui.html>
- <https://docs.maap-project.org/en/latest/system_reference_guide/jobs_maappy.html>
- <https://docs.maap-project.org/en/latest/system_reference_guide/dps_queues.html>
- <https://docs.maap-project.org/en/latest/system_reference_guide/organizations.html>
- <https://docs.maap-project.org/en/latest/system_reference_guide/personal_access_tokens.html>
- <https://docs.maap-project.org/en/latest/system_reference_guide/accessing_bucket_data.html>
- <https://docs.maap-project.org/en/latest/system_reference_guide/custom-environments.html>
- <https://docs.maap-project.org/en/latest/system_reference_guide/ade_custom_extensions/maap_libs.html>
- <https://docs.maap-project.org/en/latest/system_reference_guide/faq/ade_to_hub.html>
- <https://docs.maap-project.org/en/latest/system_reference_guide/faq/delete_algorithm_from_mas.html>
- <https://docs.maap-project.org/en/latest/troubleshooting_guides.html>
- <https://docs.maap-project.org/en/latest/troubleshooting/dps_error_codes.html>
- <https://docs.maap-project.org/en/latest/release_notes.html>
- <https://raw.githubusercontent.com/MAAP-Project/dps-unit-test/main/algorithm_config.yaml>

Non-doc sources used where the docs are silent, all marked SOURCE above:
`/srv/conda/envs/notebook/lib/python3.13/site-packages/maap/` (maap-py 4.2.0) --
`maap.py`, `AWS.py`, `Secrets.py`, `dps/DpsHelper.py`, `dps/dps_job.py`, `dps/execute.xml`,
`utils/algorithm_utils.py`, `utils/endpoints.py`.

Sections not read (judged irrelevant): Science Examples; the visualization tutorials
(titiler-pgstac, MosaicJSON, stac_ipyleaflet, lonboard, OPERA-DISP); Working with R (all
pages); Query / GEDI Cal-Val; COPC access; EDAV WCS; LPDAAC GEDI access; the
create-datasets-for-dashboard tutorial; work_with_git; ssh; kernel_resetting;
account_not_activated.
