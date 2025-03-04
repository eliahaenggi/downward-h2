import os
import json
import matplotlib.pyplot as plt

def plot_lineplots(experiment_name="hm-v1"):

    base_dir = os.path.join("data", experiment_name)

    runs_dirs = [
        os.path.join(base_dir, d) for d in os.listdir(base_dir)
        if d.startswith("runs") and os.path.isdir(os.path.join(base_dir, d))
    ]
    if not runs_dirs:
        print(f"No directory starting with 'runs' found in {base_dir}.")
        return

    all_subdirs = []
    for runs_dir in runs_dirs:
        subdirs = [os.path.join(runs_dir, subdir) for subdir in os.listdir(runs_dir) if os.path.isdir(os.path.join(runs_dir, subdir))]
        all_subdirs.extend(subdirs)

    if not all_subdirs:
        print("No runs found.")
        return

    data_by_config = {}

    for subdir in all_subdirs:
        properties_file = os.path.join(subdir, "properties")
        static_properties_file = os.path.join(subdir, "static-properties")

        if os.path.isfile(properties_file) and os.path.isfile(static_properties_file):
            with open(properties_file, "r") as f:
                properties_data = json.load(f)

            with open(static_properties_file, "r") as f:
                static_data = json.load(f)

            total_time = properties_data.get("total_time")
            config = " ".join(static_data.get("component_options", []))
            revision = static_data.get("global_revision", "unknown_revision")
            label = "astar-hmax-pi-m-compiled" if "hmax" in config else "astar-h2"
            if label not in data_by_config:
                data_by_config[label] = []
            data_by_config[label].append(total_time)

    for label in data_by_config:
        unsolved_instances = data_by_config[label].count(None)
        data_by_config[label] = sorted([time for time in data_by_config[label] if time is not None])
        data_by_config[label].extend([None] * unsolved_instances)


    output_dir = os.path.join("data", experiment_name + "-eval", "lineplots/")
    os.makedirs(output_dir, exist_ok=True)
    txt_output_path = os.path.join(output_dir, "lineplot_data.txt")

    plt.figure(figsize=(12, 8))

    with open(txt_output_path, mode='w') as txt_file:
        txt_file.write("Configuration\tTotal Time\tPercentage Solved\n")

        for label, total_times in data_by_config.items():
            solved_times = [time for time in total_times if time is not None]
            percentages = [(i + 1) / len(total_times) * 100 for i in range(len(total_times))]
            plt.plot(solved_times, percentages[:len(solved_times)], linestyle='-', label=label)

            for time, percentage in zip(solved_times, percentages[:len(solved_times)]):
                txt_file.write(f"{label}\t{time}\t{percentage}\n")

    plt.xscale('log')
    plt.title("Total Time vs Percentage of Solved Instances by Implementation", fontsize=14)
    plt.xlabel("Total Time (s)", fontsize=12)
    plt.ylabel("Percentage of Solved Instances (%)", fontsize=12)
    plt.grid(True, linestyle='-', alpha=0.7)
    plt.legend(fontsize=10, loc='best')

    output_dir = os.path.join("data", experiment_name + "-eval", "lineplots/")
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "lineplot.png")
    plt.savefig(output_path)
    plt.close()

    print(f"Lineplot saved at {output_path}")



def plot_lineplots_from_file(txt_file_path="lineplot_data.txt"):
    data_by_config = {}

    with open(txt_file_path, "r") as txt_file:
        lines = txt_file.readlines()

    for line in lines[1:]:

        parts = line.strip().split("\t")
        if len(parts) != 3:
            continue
        config, total_time, percentage_solved = parts
        total_time = float(total_time)
        percentage_solved = float(percentage_solved)

        if config not in data_by_config:
            data_by_config[config] = {'total_time': [], 'solved_instances': []}

        data_by_config[config]['total_time'].append(total_time)
        data_by_config[config]['solved_instances'].append(percentage_solved * 18.27)

    plt.figure(figsize=(12, 8))

    for config, data in data_by_config.items():
        label = config
        # dirty fix for hm and h2
        if label == "astar-h2":
            label = r"$A^* \text{ with } h^2$"
        elif label == "astar-hm":
            label = r"$A^* \text{ with } h^m$"
        plt.plot(data['total_time'], data['solved_instances'], linestyle='-', label=label)

    plt.xscale('log')
    plt.title("Total Time vs Number of Solved Instances by Implementation", fontsize=14)
    plt.xlabel("Total Time (s)", fontsize=12)
    plt.ylabel("Number of Solved Instances", fontsize=12)
    plt.grid(True, linestyle='-', alpha=0.7)
    plt.legend(fontsize=10, loc='best')

    output_dir = "lineplots/"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "lineplot_absolute_values.png")
    plt.savefig(output_path)
    plt.close()

    print(f"Lineplot saved at {output_path}")



plot_lineplots_from_file("lineplot_data.txt")
#plot_lineplots()