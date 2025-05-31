ns = [1,3]
for i in range(21):
    ns = ns+[ns[-1]*2]
ps = [1,2,8,24,48,96]

with open("results.csv", "w") as f:
    for p in ps:
        f.write(",")
        f.write(str(p))
    f.write("\n")
    for i,n in enumerate(ns):
        f.write(str(n))
        for j,p in enumerate(ps):
            filepath = "outputs/"+str(n)+"_"+str(p)+".txt"
            with open(filepath) as g:
                for line in g:
                    pass
                last_line = line
                last_line = last_line.strip()
                f.write(",")
                f.write(last_line)
        f.write("\n")