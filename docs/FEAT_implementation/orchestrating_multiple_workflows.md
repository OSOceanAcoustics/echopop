(orchestrating-multiple-workflows)=
# Orchestrate multiple workflows

Orchestration is the layer that coordinates how multiple workflows are launched, parameterized, monitored, and repeated. In practice, this usually means selecting one entrypoint pattern and then standardizing how parameters and outputs move through that pattern. The examples below are intentionally general-purpose so they can be adapted to many projects.

## Why orchestrate workflows?

When workflow logic is fragmented across various scripts or notebooks, manual execution can quickly become inconsistent. Orchestration solves this by making run behavior explicit (i.e., a run contract).

- **Reliability**: Explicitly defines the order of operations.
- **Traceability**: Captures exactly which parameters produced which outputs.
- **Scalability**: Transitions a project from "it can run on a 🥔" to "it works on a cluster".

:::{admonition} Practical framing
:class: tip
Treat orchestration as an interface design problem. The most effective setups reduce how many commands users must remember (i.e., "cognitive load") while increasing how much run context is captured automatically.
:::

(centralized-python-cli-dispatcher)=
## Centralized Python CLI dispatcher

A centralized dispatcher can be thought of as a "traffic controller". It provides a single stable CLI entrypoint that decides which downstream workflow script(s) to run. This pattern is useful when multiple workflow scripts follow similar conventions but a single stable command interface is desired for end-users. Rather than asking users to memorize each script's arguments, the dispatcher provides a single CLI surface and routes execution to the selected workflow target. Consider a case where a data workflow is modularized into three stages:

1. `workflow_ingest.py`: Data ingestion configured by `ingest.yaml`.
2. `workflow_estimation.py`: Population estimation configured by `estimate.yaml`.
3. `workflow_report.py`: Results and report generation configured by `report.yaml`.

A single dispatcher file (`workflow_dispatcher.py`) is created to coordinate separate workflow files that contain their own CLI entrypoints:

```text
project_root/
  workflows/
    workflow_ingest.py
    workflow_estimate.py
    workflow_report.py
  workflow_dispatcher.py
  configs/
    ingest.yaml
    estimate.yaml
    report.yaml
```

```{image} ../_static/orchestration_central_cli.svg
:width: 100%
:align: center
:target: ../_static/orchestration_central_cli.svg
```

When running separate workflows, the user must know the specific filename and path for every stage of the analysis. The dispatcher reduces this burden by only requiring the user to know the name of the entrypoint script and keyword for the desired stage.

::::{tab-set}
:::{tab-item} Separate workflows

Running scripts individually requires manual navigation of the `workflows/` directory.

```text
python workflows/workflow_ingest.py --config configs/ingest.yaml
python workflows/workflow_estimate.py --config configs/estimate.yaml
python workflows/workflow_report.py --config configs/report.yaml
```


:::
:::{tab-item} Dispatcher workflow

The dispatcher abstracts the directory structure, providing a single, stable command-line surface.

```powershell
python workflow_dispatcher.py --workflow ingest --config configs/ingest.yaml
python workflow_dispatcher.py --workflow estimate --config configs/estimate.yaml
python workflow_dispatcher.py --workflow report --config configs/report.yaml
```

:::
::::

:::{admonition} The "Passthrough" Benefit
:class: tip
The dispatcher can use `parse_known_args()`, meaning workflow-specific flags can still be supplied. For example:

```text
python workflow_dispatcher.py --workflow ingest --source database
```

The dispatcher ignores the `--source` flag and hands it directly to the ingestion script.
:::

### Setting it up

Only high-level control arguments should be defined within the dispatcher file, e.g., workflow target, configuration path, and other arguments. The dispatcher should also handle its own flags (e.g., which workflow to run) while passing all other "extra" flags directly to the child script. In this example, the dispatcher maintains an explicit mapping of target names to script paths.

:::{dropdown} **`workflow_dispatcher.py`**
:color: primary
:icon: code


```python
# Import dependencies
import argparse
import subprocess
import sys

# Create mapping for [workflow name ----> workflow script file]
WORKFLOW_MAP = {
    "ingest": "workflow_ingest.py",
    "estimate": "workflow_estimate.py",
    "report": "workflow_report.py",
}

# Create parser
parser = argparse.ArgumentParser(description="Central workflow dispatcher")
parser.add_argument("--workflow", required=True, choices=WORKFLOW_MAP)
parser.add_argument("--config", help="Path to config file")
args, passthrough = parser.parse_known_args()

# Build command safely
target_script = WORKFLOW_MAP[args.workflow]
cmd = [sys.executable, target_script]
if args.config:
    cmd.extend(["--config", args.config])
cmd.extend(passthrough)

result = subprocess.run(cmd, check=False)

# Ensure the dispatcher exits with the same code as the child
raise SystemExit(result.returncode)
```
:::

In production contexts, it is helpful to keep dispatcher mappings version controlled, fail fast on unknown target names, and standardize optional flags (such as verbose logging) across all child workflows.

:::{admonition} Validation checkpoint
:class: note
Before adopting this pattern broadly, test one successful run and one intentionally failing run to confirm that return codes and error messages propagate correctly through the dispatcher.
:::

### Surfacing `logging.info` from child workflow scripts

Many workflow scripts use `logging.info(...)` to announce stage boundaries, parameter choices, and summary metrics. Those messages are most useful when orchestration layers preserve them verbatim so users can see deep workflow progress at the top-level command surface.

For instance:

```python
import logging

logging.basicConfig(level=logging.INFO, format="%(message)s")

logging.info("Starting ingestion stage")
logging.info("Loading configuration from configs/ingest.yaml")
logging.info("Ingestion complete")
```

For dispatcher-style orchestration, a common pattern is to stream child process `stdout`/`stderr` directly and keep child return codes intact:

```python
import subprocess
import sys

cmd = [sys.executable, target_script, "--config", config_path]
process = subprocess.Popen(
  cmd,
  stdout=subprocess.PIPE,
  stderr=subprocess.STDOUT,
  text=True,
)

for line in process.stdout:
  print(line, end="")

raise SystemExit(process.wait())
```

With this approach, messages such as `logging.info("Starting ingestion stage")` and `logging.info("Ingestion complete")` appear in real time in the same terminal session as the orchestration command.

:::{admonition} Logging consistency tip
:class: tip
If multiple scripts are orchestrated together, keep log formatting conventions consistent (for example stage prefixes and summary blocks) so top-level logs remain readable and easy to scan.
:::

---

## Jupyter notebook

A Jupyter notebook can serve as an interactive runbook. This pattern works well when users need an interactive runbook that combines executable code, narrative notes, inline diagnostics, and run artifacts in one place. A notebook can orchestrate both Python scripts and other notebooks, which makes it useful for exploratory workflows, quality assurance runs, and tutorial-driven operations. Consider [the previous workflow example](#centralized-python-cli-dispatcher), where a single dispatcher was defined to coordinate multiple child workflow scripts. The file directory is further updated to now include `workflow_runbook.ipynb`. Consider the updated file directory:

```text
project_root/
  workflow_runbook.ipynb    <-- New orchestration runbook
  workflows/
    workflow_ingest.py
    workflow_estimate.py
    workflow_report.py
  workflow_dispatcher.py
  configs/
    ingest.yaml
    estimate.yaml
    report.yaml
  outputs/
    executed_notebooks/      <-- Papermill artifacts
```

```{image} ../_static/orchestration_notebook_hub.svg
:width: 100%
:align: center
:target: ../_static/orchestration_notebook_hub.svg
```

The orchestrator notebook acts similarly to a dispatcher script with the addition of directing the flow of data through two primary execution paths:

1. **Script execution**: The notebook calls standalone Python scripts (or a dispatcher) to perform data processing or ingestion.
2. **Child notebook execution**: The orchestrator can trigger specialized notebooks using tools like [Papermill](https://papermill.readthedocs.io/en/latest/) to generate parameterized reports or interactive analyses.

In both parts, the results are funneled into a collected outputs layer. This layer is not just a destination; it creates a feedback loop. These outputs (e.g., diagnostic plots, summary statistics) help determine if the parameters in the execution layer need adjustment. Users can return to the orchestrator to refine the configuration and re-run the stage if a result is erroneous or otherwise suspect.

### Setting up Python scripts

To trigger a workflow stage from within a notebook cell, the `subprocess` module is used. This ensures the workflow runs in its own clean process while the notebook remains in control of the sequence.

```python
# Import dependencies
import subprocess
import sys

# Trigger the dispatcher for the ingestion stage
cmd = [
    sys.executable, "workflow_dispatcher.py",
    "--workflow", "ingest",
    "--config", "configs/ingest.yaml"
]

# check=True will raise an error in the notebook if the script fails
subprocess.run(cmd, check=True)
```

### Parameterizing notebooks

Notebooks can also orchestrate *other* notebooks. Using the **Papermill** library, a notebook is treated as a function: parameters are passed in, the notebook executes, and a timestamped "output" notebook is saved as a permanent record of the run.

```python
# Import dependencies
import papermill as pm

# Generate a unique filename for the specific run
timestamp = datetime.now().strftime("%Y%m%d_%H%M")
output_path = f"outputs/executed_notebooks/ingest_run_{timestamp}.ipynb"

# Execute notebook
pm.execute_notebook(
    "input_notebook.ipynb",
    "outputs/input_notebook.executed.ipynb",
    parameters={"config_path": "configs/ingest.yaml"},
)
```

When using notebooks as orchestrators, cells **must** be designed so that repeated execution does not result in errors or redundant data processing. It is also good practice to derive relative paths from the project root to ensure scripts and configurations are located reliably. Executed notebook outputs should be written to a dedicated folder to maintain a clean source directory. Persisting these executed outputs with timestamps creates a traceable history for every analysis.

---

## Batch and shell wrappers for repeated runs

Batch and shell wrappers are lightweight scripts designed to automate the repeated execution of Python entrypoints. This avoids implementing complex looping logic within Python by leveraging the native shell environment to manage the execution sequence. The execution state is handled by the operating system, which ensures that each workflow run starts in a fresh Python process.

Consider the updated file directory where batch (`.bat`), PowerShell (`.ps1`), and shell (`.sh`) files have been added to the new folder `scripts/`:

```text
project_root/
  scripts/
    run_all_stages.bat        <-- Windows Command Prompt (CMD)
    run_all_stages.ps1        <-- Windows PowerShell
    run_all_stages.sh         <-- Linux or macOS Bash
  workflow_dispatcher.py
  workflows/
    workflow_ingest.py
    workflow_estimate.py
    workflow_report.py
  configs/
    ingest.yaml
    estimate.yaml
    report.yaml
```

```{image} ../_static/orchestration_batch_wrapper.svg
:width: 100%
:align: center
:target: ../_static/orchestration_batch_wrapper.svg
```

The orchestration logic within these scripts follows a structured execution path from initialization to result capture:

1. **Wrapper loop logic**: The script (e.g., `.bat`, `.ps1`, or `.sh`) serves as the primary controller. It defines a workflow list or array that iterates through specific targets, such as different analyses or survey-specific processing files. For each target, the wrapper prepares the environment, sets necessary variables, and constructs the execution command.

2. **Python entrypoint**: The wrapper invokes the Python dispatcher or a specific script as a child process. By passing a pre-run configuration (assuming `--config` is still used as a file pointer) and any captured passthrough arguments, the wrapper ensures that each run is independent and reproducible.

3. **Result capture and logging**: As each process completes, the wrapper captures the exit code (e.g., `$LASTEXITCODE` or `errorlevel`). This determines whether the run was a success or a failure. These results are funneled into a summary or dedicated log file, providing an audit trail for the entire batch.

::::{tab-set}
:::{tab-item} `.bat`

Windows command prompt loop.

```bat
@echo off
setlocal enabledelayedexpansion

:: Iterate through the defined workflow names
for %%W in (ingest estimate report) do (
    echo [INFO] Executing workflow: %%W ...

    :: Call the dispatcher from the parent directory
    python ..\workflow_dispatcher.py --workflow %%W --config ..\configs\%%W.yaml

    :: Check the exit code of the Python process
    if errorlevel 1 (
        echo [ERROR] Workflow %%W failed.
    )
)
```

This batch file can then be run in the Windows Command Prompt (`cmd.exe`) via:

```bat
cd scripts
run_all_stages.bat
```

:::
:::{tab-item} `.ps1`

PowerShell loop.

```powershell
$workflows = @("ingest", "estimate", "report")

foreach ($w in $workflows) {
    Write-Host ">>> Starting $w stage..." -ForegroundColor Cyan

    python ../workflow_dispatcher.py --workflow $w --config "../configs/$w.yaml"

    if ($LASTEXITCODE -ne 0) {
        Write-Warning "Execution failed for $w. Check logs for details."
    }
}
```

This file can be run in PowerShell via:

```powershell
cd scripts
.\run_all_stages.ps1
```

Note that PowerShell has stricter security policies than CMD. By default, it may prevent the execution of local scripts unless the "Execution Policy" is adjusted *for the current session*. This requires updating the command to:

```powershell
cd scripts

# Run with temporary policy bypass if restricted
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass

# Execute the script
.\run_all_stages.ps1
```

:::
:::{tab-item} `.sh`

Linux/macOS shell loop.

```bash
#!/bin/bash

# Define array of stages
workflows=(ingest estimate report)

for w in "${workflows[@]}"; do
    echo "--- Running $w ---"

    # Run python and echo failure if the exit code is non-zero
    python3 ../workflow_dispatcher.py --workflow "$w" --config "../configs/${w}.yaml" || \
    echo "[CRITICAL] $w stage failed"
done
```

In Unix-like systems, a script must be granted "execute" permissions before it can be run. The resulting command becomes:

```bash
cd scripts

# Grant execution permission (only needed once)
chmod +x run_many_workflows.sh

# Execute the script
./run_many_workflows.sh
```

:::
::::

Batch wrappers are a practical option when you need repeated execution over a fixed list of targets and want minimal overhead. Instead of building a full orchestration service, users define a loop in a shell wrapper and call one Python entrypoint for each target run.

---

## Task runner entrypoints

Task runners provide a command vocabulary for a project. Instead of requiring users to remember all flags, directory paths, or specific Python interpreters, task runners allow for the definition of short, memorable aliases (e.g., `do ingest`, `kriging population`). This pattern is especially valuable in collaborative environments. It serves as a form of executable documentation, where the standard operating procedures for the project are encoded in a single, version-controlled file (e.g., `Makefile`, `justfile`).

Consider the updated file directory where batch (`.bat`), PowerShell (`.ps1`), and shell (`.sh`) files have been added to the new folder `scripts/`:

```text
project_root/
  Makefile or justfile        <-- Command definitions
  workflows/
    workflow_ingest.py
    workflow_estimate.py
    workflow_report.py
  workflow_dispatcher.py
  configs/
    ingest.yaml
    estimate.yaml
    report.yaml
```

```{image} ../_static/orchestration_task_runner.svg
:width: 100%
:align: center
:target: ../_static/orchestration_task_runner.svg
```

### Using task runners

Task runners provide named commands for repeat orchestration steps. This is useful for teams that want a clean command vocabulary (e.g., `run-ingest`, `run-estimate`, `run-report`, `run-all`) without asking users to remember long command lines. The task runner pattern acts as a user interface layer for the command line. It translates simple, human-readable instructions into the technical, multi-argument commands required by the underlying Python workflows.

#### Task runner

This is the top-level interface that represents the user experience. In this layer, the complexity of file paths and flag names is hidden. A user only needs to interact with a small set of standardized commands, e.g.,

- `make run-ingest`
- `make run-estimate`
- `make run-report`
- `make run-all`

#### Standard command layer

When a task is triggered, the runner invokes the standard command layer. This is where the run contract is enforced. The runner takes the short alias and expands it into the full execution string. This layer ensures that every time a task is run, it uses the exact same dispatcher and the exact same configuration pointer, eliminating manual entry errors.

::::{tab-set}
:::{tab-item} Named task

```text
make run-ingest
```

:::
:::{tab-item} Dispatcher

```text
python workflow_dispatcher.py --workflow ingest --config configs/ingest.yaml
```
:::
::::

#### Workflow scripts

The dispatcher routes the command to the specific workflow scripts. At this stage, the domain-specific logic is executed. Because this layer is separated from the interface layer, the actual Python code can be updated or refactored without changing the commands the user types into the terminal.

#### Output

The execution produces measurable results: processed data artifacts, execution logs, and summary reports.

### Minimal `Makefile` example

The simplest implementation of a `Makefile` acts as a direct alias for the dispatcher commands. In this form, each "target" (i.e., the word preceding each colon) represents a named task that executes a specific command.

```make
run-ingest:
	python workflow_dispatcher.py --workflow ingest --config configs/ingest.yaml

run-estimate:
	python workflow_dispatcher.py --workflow estimate --config configs/estimate.yaml

run-report:
    python workflow_dispatcher.py --workflow report --config configs/report.yaml
```

### Implementing chain dependencies

One of the primary advantages of a task runner is the ability to define dependencies. In scientific workflows, certain stages should only execute if the previous stage succeeded. By listing other targets on the same line as a new target, a "chain" is created. This adds a new task called `run-all` to the `Makefile`:

```make
run-ingest:
	python workflow_dispatcher.py --workflow ingest --config configs/ingest.yaml

run-estimate:
	python workflow_dispatcher.py --workflow estimate --config configs/estimate.yaml

run-report:
    python workflow_dispatcher.py --workflow report --config configs/report.yaml

run-all: run-ingest run-estimate run-report
```

When `make run-all` is triggered, the runner executes the tasks from left to right. A critical feature of this logic is the fail-fast mechanism: if `run-ingest` returns a non-zero exit code (indicating a crash or error), the runner halts the entire process. This prevents the "estimate" stage from running on corrupted or missing data.

### Advanced task runner logic

As a research project grows, the task runner is often expanded to include shortcuts (variables) and cleanup steps. This ensures that the entire team is using the exact same version of the code and the same configuration files, even if individuals are working on different computers or updated software. By centralizing these commands, the "recipe" for the analysis stays consistent, preventing the accidental use of the wrong data versions or outdated settings.

To maintain this consistency as the project gains complexity, several specific features of the task runner are utilized:

- **Variable management**: Definitions like `PYTHON = python` at the top of the file act as global constants. If the project migrates to a newer version of Python or moves the configuration folder, the update is made in one line rather than across every individual task.

- **Phony targets**: The `.PHONY` instruction tells the computer that these names are commands to be executed, not actual files on the hard drive. This prevents the runner from getting confused if a folder named `clean` or `report` happens to exist in the project directory.

- **Cleanup tasks**: A `clean` task uses system commands to wipe old logs or temporary files. This ensures the results folder is reset to a "clean slate" before a fresh analysis begins, preventing data from previous runs from contaminating new results.

:::{dropdown} **`Makefile`**
:color: primary
:icon: code

```make
# Variables for consistency
PYTHON = python
DISPATCH = workflow_dispatcher.py
CONF_DIR = configs

.PHONY: run-ingest run-estimate run-report run-all clean

# Individual stage tasks
run-ingest:
	$(PYTHON) $(DISPATCH) --workflow ingest --config $(CONF_DIR)/ingest.yaml

run-estimate:
	$(PYTHON) $(DISPATCH) --workflow estimate --config $(CONF_DIR)/estimate.yaml

run-report:
	$(PYTHON) $(DISPATCH) --workflow report --config $(CONF_DIR)/report.yaml

# Chained execution
run-all: run-ingest run-estimate run-report
	@echo "Full pipeline successfully executed."

# Utility task
clean:
	rm -rf outputs/logs/*.log
	rm -rf __pycache__
```
:::

::::::{admonition} Cross-platform note for `make`
:class: tip
`make` is typically available by default on Linux and macOS. On Windows, it is usable but usually requires installing GNU Make (for example through MSYS2, Git Bash/MinGW, Cygwin, WSL, or Chocolatey). Other common tool choices include [`just`](https://just.systems/man/en/), [`Invoke`](https://www.pyinvoke.org/), and [`Nox`](https://nox.thea.codes/en/stable/), depending on team preference and platform constraints. The above `Makefile` can be readily translated to any of these.

:::::{tab-set}
::::{tab-item} `just`

The `just` runner is purpose-built for command shortcuts and works natively across Windows, macOS, and Linux without the strict formatting requirements of `make`. It uses a syntax very similar to `make` but handles variables and cross-platform pathing more cleanly.

:::{dropdown} Click to expand `justfile` code
:color: light

```make
# Variables for consistency
python := "python"
dispatch := "workflow_dispatcher.py"
conf_dir := "configs"

# Individual stage tasks
run-ingest:
    {{python}} {{dispatch}} --workflow ingest --config {{conf_dir}}/ingest.yaml

run-estimate:
    {{python}} {{dispatch}} --workflow estimate --config {{conf_dir}}/estimate.yaml

run-report:
    {{python}} {{dispatch}} --workflow report --config {{conf_dir}}/report.yaml

# Chained execution
run-all: run-ingest run-estimate run-report
    @echo "Full pipeline successfully executed."

# Utility task
clean:
    rm -rf outputs/logs/*.log
    rm -rf __pycache__
```
:::
::::
::::{tab-item} `Invoke`

`Invoke` is a Python-native task runner. It is highly effective for teams that prefer to stay entirely within the Python ecosystem, as it allows for the use of standard Python logic to handle file paths and execution.

:::{dropdown} Click to expand tasks.py code
:color: light
```python
from invoke import task

# Variables for consistency
PYTHON = "python"
DISPATCH = "workflow_dispatcher.py"
CONF_DIR = "configs"

# Individual stage tasks
@task
def run_ingest(c):
    c.run(f"{PYTHON} {DISPATCH} --workflow ingest --config {CONF_DIR}/ingest.yaml")

@task
def run_estimate(c):
    c.run(f"{PYTHON} {DISPATCH} --workflow estimate --config {CONF_DIR}/estimate.yaml")

@task
def run_report(c):
    c.run(f"{PYTHON} {DISPATCH} --workflow report --config {CONF_DIR}/report.yaml")

# Chained execution
@task(pre=[run_ingest, run_estimate, run_report])
def run_all(c):
    print("Full pipeline successfully executed.")

# Utility task
@task
def clean(c):
    c.run("rm -rf outputs/logs/*.log")
    c.run("rm -rf __pycache__")
```
:::
::::
::::{tab-item} `Nox`
`Nox` focuses on reproducibility by automatically creating isolated environments for every task. This ensures the analysis is not dependent on specific packages installed on one person's laptop.

:::{dropdown} Click to expand noxfile.py code
:color: light
```python
import nox

# Variables for consistency
PYTHON = "python"
DISPATCH = "workflow_dispatcher.py"
CONF_DIR = "configs"

# Individual stage tasks
@nox.session
def run_ingest(session):
    session.run(PYTHON, DISPATCH, "--workflow", "ingest", "--config", f"{CONF_DIR}/ingest.yaml")

@nox.session
def run_estimate(session):
    session.run(PYTHON, DISPATCH, "--workflow", "estimate", "--config", f"{CONF_DIR}/estimate.yaml")

@nox.session
def run_report(session):
    session.run(PYTHON, DISPATCH, "--workflow", "report", "--config", f"{CONF_DIR}/report.yaml")

# Chained execution
@nox.session
def run_all(session):
    session.notify("run_ingest")
    session.notify("run_estimate")
    session.notify("run_report")

# Utility task
@nox.session
def clean(session):
    session.run("rm", "-rf", "outputs/logs/*.log", external=True)
    session.run("rm", "-rf", "__pycache__", external=True)
```
:::
::::
:::::
::::::

---

## Choosing a pattern

In many real-world deployments, these patterns are combined rather than used in isolation. A common arrangement is a **nested hierarchy**: a scheduler or batch script triggers a task runner, which calls the centralized dispatcher, which finally executes the individual workflow scripts. Selecting a pattern should be based on who is running the workflows, the required frequency of runs, and the level of automation needed.

### Selecting an approach

The following decision matrix summarizes the trade-offs between the orchestration patterns to assist in selecting an appropriate entrypoint. The above examples are also not an exhaustive representation of all possible patterns.

| Pattern | Best for... | User Experience | Primary Strength |
| :---- | :------- | :---- | :---- |
| **Dispatcher** | Power users & developers | CLI-driven; requires knowing flags. | Single "source of truth" for script routing. |
| **Notebook** | Exploration & reporting | Interactive; visual and narrative. | Documentation and results live side-by-side. |
| **Batch wrapper** | High-volume repetition | Hands-off; automated loops. | Leverages the OS to run many years/surveys. |
| **Task runner** | Collaborative teams | Simple aliases (e.g., `make run-all`). | Executable documentation; manages dependencies. |

The choice of pattern is often driven by specific operational goals of the research:

- **Standardization**: If the primary goal is to ensure every user passes the same flags and configuration files to a set of fragmented scripts, the **Dispatcher** pattern is the most direct solution.
- **Ease of use**: If the project is shared among a team with varying levels of command-line experience, a **Task runner** provides a simplified command vocabulary that hides technical complexity.
- **Reporting and transparency**: When the final output must be an interactive or visual document that explains the logic of the analysis, a **Notebook**-based orchestrator may be preferred.
- **Scale**: When an analysis must be repeated across dozens of survey years or geographic regions without manual intervention, a **Batch wrapper** provides the necessary looping logic.

:::{admonition} Setting users up for success
:class: tip
Regardless of the pattern selected, the objective of a workflow orchestrator should be to reduce the cognitive load on end-users while enhancing scientific traceability. A robust setup ensures that every result in the `outputs/` folder (or elsewhere) can be traced back to the exact command and configuration file that produced it.
:::
