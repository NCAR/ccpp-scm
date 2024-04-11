import argparse
import itertools
import os
import subprocess
import sys
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# TODO:
#  - [ ] add user argument to cmake "$@" in configure{32,64}
#  - [ ] add run function
#  - [ ] run through all cases compiled
suites = ['SCM_RAP']
cases = ['twpice']

global verbose
verbose = False

def print_cmd(command):
    print(f"$ {' '.join(command)}")

def run_cmd(command):
    print_cmd(command)
    try:
        # Run the cmake command
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        # Handle errors if cmake command fails
        print(f"Error running {e}")
        sys.exit()


# Function to run configuration step
def configure32():
    print("Running single precision configuration step")
    command = ["cmake", "src/", "-B", "build32", "-D32BIT=1"]
    run_cmd(command)


# Function to run configuration step
def configure64():
    print("Running double precision configuration step")
    command = ["cmake", "src/", "-B", "build64"]
    run_cmd(command)


# Function to run configuration step
def configure():
    print("Running configuration steps")
    configure32()
    configure64()


# Function to run build step
def build32():
    print("Running single precision build step")
    command = ["make", "-C", "build32", "-j4"]
    run_cmd(command)

# Function to run build step
def build64():
    print("Running double precision build step")
    command = ["make", "-C", "build64", "-j4"]
    run_cmd(command)


# Function to run build step
def build():
    print("Running make steps")
    build32()
    build64()


# Post process single vs. double precision output
def post():
    print("Post processing single vs. double precision output")
    postprocess()


def run32():
    print("Running single precision executable")
    combinations = list(itertools.product(cases, suites))
    for (case, suite) in combinations:
        command = ["build32/run_scm.py", "-c", case,
                   "-s", suite,
                   "--bin_dir", os.environ["PWD"]+"/build32",
                   "--run_dir", os.environ["PWD"]+"/output32"]
        run_cmd(command)

def run64():
    print("Running double precision executable")
    combinations = list(itertools.product(cases, suites))
    for (case, suite) in combinations:
        command = ["build64/run_scm.py", "-c", case,
                   "-s", suite,
                   "--bin_dir", os.environ["PWD"]+"/build64",
                   "--run_dir", os.environ["PWD"]+"/output64"]
        run_cmd(command)

# Function to run executable
def run():
    print("Running executables")
    run32()
    run64()

def check_for_vertical_dimension(da):
    vert_dim = None
    if 'vert_dim_layer' in da.dims:
        vert_dim = True
    return vert_dim


def get_time_dimension(da):
    # Initialize time variable to None
    time = None

    # Iterate over dimensions of the Dataset
    for dim in da.dims:
        if 'time' in dim.lower():
            time = da[dim].name
            break

    # Check if time variable was found
    if time is None:
        print(f"No 'time' dimension matched in {da.name}")
    # else:
    #     print("Found 'time' dimension in the dataset.")
    return time

def check_dimensions(da1):
    found_time = True
    found_vert_layer = True
    time_dim = get_time_dimension(da1)
    vert_dim = check_for_vertical_dimension(da1)
    if time_dim == None:
        found_time = False
        if verbose:
            print(f"No time dimension found in variable {da1.name}")
        return found_time
    if vert_dim == None:
        found_vert_layer = False
        if verbose:
            print(f"No vertical layer found in variable {da1.name}")
        return found_vert_layer
    return time_dim

def calculate_mse(da1, da2, time_dim):
    print("TIME DIM", time_dim)
    """Calculate the mean square error (MSE) between two DataArrays."""
    mse = ((da1 - da2) ** 2)#.mean(dim=time_dim)
    return mse


def save_mse_plots(case_dir, mse_dict):
    """Save plots of mean square error (MSE) to a directory."""
    diff_dir = case_dir + '_diff'
    os.makedirs(diff_dir, exist_ok=True)

    for var_name, mse in mse_dict.items():
        plt.figure()
        mse.plot()
        plt.title(f'Mean Square Error for {var_name}')
        plt.xlabel('Column Level')
        plt.ylabel('MSE')
        plt.savefig(os.path.join(diff_dir, f'{var_name}_mse.png'))
        plt.close()

def save_re_plots(case_dir, re_dict):
    """Save plots of relative error to a directory."""
    diff_dir = case_dir + '_diff'
    os.makedirs(diff_dir, exist_ok=True)

    for var_name, re in re_dict.items():
        plt.figure()
        re.plot()
        plt.title(f'Relative Error for {var_name}')
        plt.xlabel('Column Level')
        plt.ylabel('Relative Error')
        plt.savefig(os.path.join(diff_dir, f'{var_name}_re.png'))
        plt.close()


def save_percentage_plots(case_dir, pd_dict):
    """Save plots of percentage difference to a directory."""
    diff_dir = case_dir + '_diff'
    os.makedirs(diff_dir, exist_ok=True)

    for var_name, pd in pd_dict.items():
        plt.figure()
        pd.plot()
        plt.title(f'Percentage Difference for {var_name}')
        plt.xlabel('Column Level')
        plt.ylabel('Percentage Diff (%)')
        plt.savefig(os.path.join(diff_dir, f'{var_name}_pd.png'))
        plt.close()

# During post-processing step calculate relative the percentage difference
def calculate_relative_error_and_percentage(da32, da64, time_dim):
    """Calculate the relative error between two DataArrays."""
    found_time = True
    time_dim = get_time_dimension(da32)
    if time_dim == None:
        found_time = False
        print(f"No time dimension found in variable {da32.name}")
        return None, found_time

    relative_error = (np.abs(da64 - da32) / da64).mean(dim=time_dim)
    percentage_dif = relative_error * 100
    return relative_error, percentage_dif


# Post processing function
def postprocess():
    print('Post-processing')
    analysis_dir = 'output_diff'
    analysis_file = 'variable_report.txt'
    output_f = 'output.nc'
    output32_dir = 'output32'
    output64_dir = 'output64'
    output32_dirs = sorted([os.path.basename(dirpath) for dirpath, _, _ in os.walk(output32_dir) if os.path.basename(dirpath).startswith('output_')])
    output64_dirs = sorted([os.path.basename(dirpath) for dirpath, _, _ in os.walk(output64_dir) if os.path.basename(dirpath).startswith('output_')])

    vars_equal = []
    vars_diff = []

    # Use set intersection to find directories that exist in both output32 and output64
    case_dirs = set(output32_dirs).intersection(output64_dirs)
    for case_dir in case_dirs:
        f32 = os.path.join(output32_dir, case_dir, output_f)
        f64 = os.path.join(output64_dir, case_dir, output_f)
        # Open f32 and f64 using xarray
        ds32 = xr.open_dataset(f32)
        ds64 = xr.open_dataset(f64)

        # Compare variables and compute mean square error (MSE) for variables that are different
        mse_dict = {}
        re_dict = {}
        percentage_diff_dict = {}
        for var_name in ds32.data_vars:
            if var_name in ds64.data_vars:
                if np.any(ds32[var_name] != ds64[var_name]):
                    time_dim = check_dimensions(ds32[var_name])
                    if time_dim == False:
                        continue
                    mse = calculate_mse(ds32[var_name], ds64[var_name], time_dim)
                    relative_error, percentage_diff = \
                        calculate_relative_error_and_percentage(ds32[var_name],
                                                                ds64[var_name],
                                                                time_dim)
                    mse_dict[var_name] = mse
                    re_dict[var_name] = relative_error
                    percentage_diff_dict[var_name] = percentage_diff
                    vars_diff.append(var_name)
                    diff_f = os.path.join(analysis_dir,case_dir)
                    print(f"Variable {var_name} differs in output files, writing diff to {diff_f}")
                    # Save plots of mean square error (MSE) to a directory
                    save_mse_plots(diff_f, mse_dict)
                    save_re_plots(diff_f, re_dict)
                    save_percentage_plots(diff_f, percentage_diff_dict)
                else:
                    vars_equal.append(var_name)
                    print(f"Variable {var_name} is the same in both output files")

    # Write the vars_diff list to a text file
    with open(os.path.join(analysis_dir, analysis_file), 'w') as f:
        f.write("Different variables:\n")
        for var in vars_diff:
            f.write(f" - {var}\n")
        f.write("Equal variables:\n")
        for var in vars_equal:
            f.write(f" - {var}\n")


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Script to configure, build, and run a program")
    parser.add_argument("--configure", action="store_true", help="Configure CMake for single and double precision")
    parser.add_argument("--configure32", action="store_true", help="Configure CMake for single precision")
    parser.add_argument("--configure64", action="store_true", help="Configure CMake for double precision")
    parser.add_argument("--build", action="store_true", help="Build single and double precision code")
    parser.add_argument("--build32", action="store_true", help="Build single precision code")
    parser.add_argument("--build64", action="store_true", help="Build double precision code")
    parser.add_argument("--post", action="store_true", help="Postprocess the output")
    parser.add_argument("--run", action="store_true", help="Run cases in single and double precision")
    parser.add_argument("--run32", action="store_true", help="Run cases in single precision")
    parser.add_argument("--run64", action="store_true", help="Run cases in double precision")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    args = parser.parse_args()

    # If --help is provided, print help message
    if len(sys.argv) == 1 or "--help" in sys.argv:
        parser.print_help()
        sys.exit()

    verbose = args.verbose

    # Execute requested step(s)
    if args.configure:
        configure()
    elif args.configure32:
        configure32()
    elif args.configure64:
        configure64()
    elif args.build32:
        build32()
    elif args.build64:
        build64()
    elif args.build:
        build()
    elif args.run:
        run()
    elif args.run32:
        run32()
    elif args.run64:
        run64()
    elif args.post:
        post()
    elif args.run:
        run()
    else:
        configure()
        build()
        run()


if __name__ == "__main__":
    main()
