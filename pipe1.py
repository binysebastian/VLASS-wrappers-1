import subprocess
import time
import re

def run_command(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return stdout.decode('utf-8'), stderr.decode('utf-8'), process.returncode

def is_pipeline_running():
    stdout, stderr, exit_code = run_command("python3 catenator.py monitor")
    if exit_code == 0:
        return "* * * PIPELINE IS RUNNING * * *" in stdout
    else:
        print(f"Error checking pipeline status: {stderr}")
        return False

def wait_for_pipeline():
    while is_pipeline_running():
        print("Pipeline is running... waiting for 2 minutes before checking again.")
        time.sleep(120)  # Wait for 2 minutes

commands = [
    "python3 catenator.py download",
    "python3 catenator.py process",
    "python3 catenator.py vlad",
    "python3 catenator.py diagnostics",
    "python3 catenator.py qa"
]

for command in commands:
    print(f"Running command: {command}")
    stdout, stderr, exit_code = run_command(command)
    if exit_code != 0:
        print(f"Error running command '{command}': {stderr}")
        break
    print(stdout)
    wait_for_pipeline()

print("Automation script completed.")
