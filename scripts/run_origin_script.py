"""
This script executes another Python script, `get_physical_origin_values.py`,
using a subprocess and captures its output.

The script performs the following steps:
1.  Constructs the full path to the `get_physical_origin_values.py` script,
    assuming it is located in the same directory.
2.  Uses `subprocess.run` to execute the target script.
3.  Captures and prints the standard output and standard error of the executed script.
4.  Includes error handling to catch and report any issues during script execution.
"""
import subprocess
import os

# The script is assumed to be in the 'scripts' directory.
# The target script is in the same directory.
script_dir = os.path.dirname(os.path.abspath(__file__))
script_name = "get_physical_origin_values.py"
full_script_path = os.path.join(script_dir, script_name)

try:
    result = subprocess.run(["python", full_script_path], capture_output=True, text=True, check=True)
    print(result.stdout)
    print(result.stderr)
except subprocess.CalledProcessError as e:
    print(f"Error executing script: {e}")
    print(f"Stdout: {e.stdout}")
    print(f"Stderr: {e.stderr}")
