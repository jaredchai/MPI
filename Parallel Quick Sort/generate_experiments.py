import random

def generate_random_integers(n):
    return [random.randint(-n, n) for _ in range(n)]

def save_random_integers_to_file(filepath, integers):
    with open(filepath, "w") as f:
        f.write(f"{len(integers)}\n")
        f.write(" ".join(str(i) for i in integers))

ns = [1,3]
for i in range(21):
    ns = ns+[ns[-1]*2]
ps = [1,2,8,24,48,96]
PBS = "#PBS -N pqsort\n#PBS -A tchai8\n#PBS -l nodes=4:ppn=24\n#PBS -l pmem=1gb\n#PBS -l walltime=30:00\n#PBS -q coc-ice-multi\n#PBS -j oe\n#PBS -o pqsort.outfile\n\ncd $PBS_O_WORKDIR\n\n"
for n in ns:
    random_integers = generate_random_integers(n)
    file = str(n)+".txt"
    filepath = "inputs/"+file
    save_random_integers_to_file(filepath, random_integers)
    PBS = PBS + "echo \"n = "+str(n)+":\"\n"
    for p in ps:
        PBS = PBS + "echo \"p = "+str(p)+":\"\n"
        PBS = PBS + "mpiexec -np "+str(p)+" ./pqsort "+filepath+" outputs/"+str(n)+"_"+str(p)+".txt\n"
    PBS = PBS+"\n"
with open("pqsort.pbs", "w") as f:
    f.write(PBS)