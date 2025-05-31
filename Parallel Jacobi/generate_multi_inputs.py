import subprocess

def main():
    n_values = list(map(int, input("Enter the values for n (space-separated): ").split()))

    for n in n_values:
        # Prepare the suffix for the input files
        suffix = f"{n}"
        # Call the generate_inputs.py script with the n value and suffix
        subprocess.run(["python", "generate_input.py", str(n), suffix])

    print("Input files generated.")

if __name__ == "__main__":
    main()