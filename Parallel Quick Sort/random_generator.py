import random
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate a file with a list of random integers.")
    parser.add_argument("n", type=int, help="Number of integers in the list.")
    parser.add_argument("filepath", type=str, help="File path.")
    return parser.parse_args()

def generate_random_integers(n):
    return [random.randint(-100, 100) for _ in range(n)]

def save_random_integers_to_file(filepath, integers):
    with open(filepath, "w") as f:
        f.write(f"{len(integers)}\n")
        f.write(" ".join(str(i) for i in integers))

def main():
    args = parse_arguments()
    random_integers = generate_random_integers(args.n)
    save_random_integers_to_file(args.filepath, random_integers)

if __name__ == "__main__":
    main()