import subprocess
import os
from pathlib import Path
import argparse

problem_files = {
    "gripper/domain.pddl": [
        "gripper/prob01.pddl"
    ],
    "npuzzle/n-puzzle-typed.pddl": [
        "npuzzle/p2.pddl",
        "npuzzle/p3_1.pddl",
        "npuzzle/p3_2.pddl",
        "npuzzle/p3_3.pddl",
        "npuzzle/p3_4.pddl",
        "npuzzle/p4.pddl",
    ]
}

python_executable = "python3" if os.name == "posix" else "python"

path_to_fast_downward = "../fast-downward.py"

# Default search strategy and time limit
default_search_strategy = 'astar(hm())'
default_time_limit = 60
default_log_file = "output.txt"

def run_fast_downward(domain, problem, search_strategy, time_limit):
    command = [
        python_executable,
        path_to_fast_downward,
        domain,
        problem,
        "--search",
        search_strategy
    ]
    
    if not Path(domain).exists() or not Path(problem).exists():
        return None, f"Error: File not found - Domain: {domain}, Problem: {problem}"

    try:
        result = subprocess.run(
            command, 
            check=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            timeout=time_limit
        )
        output = result.stdout.decode('utf-8')
        return output, None
    except subprocess.TimeoutExpired:
        return None, f"Timeout expired after {time_limit} seconds for Domain: {domain}, Problem: {problem}"
    except subprocess.CalledProcessError as e:
        return None, e.stderr.decode('utf-8')

def parse_output(output):
    lines = output.splitlines()
    relevant_info = {}
    for line in lines:
        if "Plan length" in line:
            relevant_info["Plan length"] = line.split(":")[1].strip()
        elif "Plan cost" in line:
            relevant_info["Plan cost"] = line.split(":")[1].strip()
        elif "Expanded " in line and "state(s)" in line and "until last jump" not in line:
            relevant_info["Expanded states"] = line.split("Expanded")[1].split("state")[0].strip()
        elif "Search time" in line:
            relevant_info["Search time"] = line.split(":")[1].strip()
        elif "Total time" in line:
            relevant_info["Total time"] = line.split(":")[1].strip()
    return relevant_info


def main():
    parser = argparse.ArgumentParser(description="Run Fast Downward with a specified search strategy and time limit.")
    parser.add_argument(
        "--search",
        type=str,
        default=default_search_strategy,
        help="Specify the search strategy (default: 'astar(hm())')"
    )
    parser.add_argument(
        "--tlimit",
        type=int,
        default=default_time_limit,
        help="Specify the time limit for each problem in seconds (default: 60)"
    )
    parser.add_argument(
        "--output",
        type=str,
        default=default_log_file,
        help="Specify the log file name (default: 'output.txt')"
    )
    args = parser.parse_args()
    search_strategy = args.search
    time_limit = args.tlimit
    log_file = args.output

    if not Path(path_to_fast_downward).exists():
        print("Error: 'fast-downward.py' not found.")
        return

    with open(log_file, "w") as log:
        log.write("Fast Downward Execution Log\n")
        log.write("=" * 40 + "\n")

        for domain in problem_files.keys():
            for problem in problem_files[domain]:
                log.write(f"Running with Domain: {domain}, Problem: {problem}, Search Strategy: {search_strategy}, Time Limit: {time_limit} seconds\n")
                print(f"Running with Domain: {domain}, Problem: {problem}, Search Strategy: {search_strategy}, Time Limit: {time_limit} seconds")

                output, error = run_fast_downward(domain, problem, search_strategy, time_limit)
                
                if error:
                    log.write("Execution failed with error:\n")
                    log.write(error + "\n")
                    print(f"Execution failed with error: {error}")
                else:
                    info = parse_output(output)
                    for key, value in info.items():
                        log.write(f"{key}: {value}\n")
                        print(f"{key}: {value}")
                log.write("-" * 40 + "\n")

if __name__ == "__main__":
    main()
