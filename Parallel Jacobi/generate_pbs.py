def generate_pbs_script(n_values, p_values):
    content = f"""#PBS -N pjacobi
#PBS -A hren62
#PBS -l nodes=4:ppn=24
#PBS -l pmem=1gb
#PBS -l walltime=30:00
#PBS -q coc-ice-multi
#PBS -j oe
#PBS -o pjacobi.outfile

cd $PBS_O_WORKDIR
"""

    for j in n_values:
        content+= f"""

echo "n = {j}:" """
        for i in p_values:
            content += f"""
echo "p = {i}:"
mpiexec -np {i} ./pjacobi input/input_mat_{j}.txt input/input_vec_{j}.txt output/out_{j}_{i}.txt"""

    return content


def main():
    n_values = list(map(int, input("Enter the values for n (space-separated): ").split()))
    p_values = list(map(int, input("Enter the values for p (space-separated) (should be perfect square): ").split()))

    script_content = generate_pbs_script(n_values, p_values)

    with open("pjacobi.pbs", "w") as f:
        f.write(script_content)

    print("PBS script generated and saved as 'pjacobi.pbs'.")


if __name__ == "__main__":
    main()
