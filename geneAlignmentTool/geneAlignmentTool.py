import heapq

class searchAgent:
    def __init__(self, gnm, gene, k=9):
        self.gnm = gnm
        self.gene = gene
        self.k = k
        self.matchScore = 0
        self.mismatchScore = 2
        self.gapScore = 3

    def getSeeds(self): # getting seeds of 9mers using a sliding window
        seed = self.gene[:self.k]
        seeds = []
        for i in range(len(self.gnm) - self.k):
            if self.gnm[i:i+self.k] == seed:
                seeds.append(i)
        return seeds

    def heuristic(self, idx): 
        return len(self.gene) - idx

    def aStarSearch(self, sidx): # heuristic search function with backtracking
        start = (sidx, 0)
        frontier = [(self.heuristic(0), 0, start)]
        pth = {start: None}
        costDict = {start: 0}

        bestState = None
        bestCost = float('inf')

        while frontier:
            f, g, cur = heapq.heappop(frontier)
            gi, gj = cur

            if gj == len(self.gene):
                if g < bestCost:
                    bestCost = g
                    bestState = cur
                break

            moves = []
            if gi < len(self.gnm) and gj < len(self.gene):
                if self.gnm[gi] == self.gene[gj]:
                    cost = self.matchScore
                else:
                    cost = self.mismatchScore
                moves.append((gi+1, gj+1, cost))
            if gj < len(self.gene):
                moves.append((gi, gj+1, self.gapScore))
            if gi < sidx + len(self.gene) + 10:
                moves.append((gi+1, gj, self.gapScore))

            for ngi, ngj, moveCost in moves:
                nxt = (ngi, ngj)
                newCost = costDict[cur] + moveCost
                if nxt not in costDict or newCost < costDict[nxt]:
                    costDict[nxt] = newCost
                    heapq.heappush(frontier, (newCost + self.heuristic(ngj), newCost, nxt))
                    pth[nxt] = cur

        if bestState is None:
            return ("No alignment", ""), float('inf')

        gnmAln = []
        geneAln = []
        cur = bestState
        while pth[cur] is not None:
            prev = pth[cur]
            if cur[0] > prev[0] and cur[1] > prev[1]:
                gnmAln.append(self.gnm[prev[0]])
                geneAln.append(self.gene[prev[1]])
            elif cur[0] > prev[0]:
                gnmAln.append(self.gnm[prev[0]])
                geneAln.append("-")
            else:
                gnmAln.append("-")
                geneAln.append(self.gene[prev[1]])
            cur = prev

        return ("".join(reversed(gnmAln)), "".join(reversed(geneAln))), bestCost

def readSeq(path):
    seq = ""
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            seq += line
    return seq

#alter the input path to enter other files
gnmPath = "GCA_000005845.2_ASM584v2_genomic.fna"
genePath = "gene.fna"

gnm = readSeq(gnmPath)
gene = readSeq(genePath)

# genome and gene sequences can be directly entered here
#gnm = "AAAAAACCCCGGGTTTTT"
#gene   = "AAAAAACCCGGGTTTTT" 

agent = searchAgent(gnm, gene)
seeds = agent.getSeeds()
print(f"No. of seeds found = {len(seeds)}")


res = []
for s in seeds:
    (a1, a2), score = agent.aStarSearch(s)
    if a1 == "No alignment":
        continue
    res.append((score, s, a1, a2))
    print(f"Processed seed {s} -> score {score}")

with open("allAlignments.txt", "w") as f:
    for score, pos, a1, a2 in res:
        f.write(f"Seed at {pos}\n")
        f.write(f"Score = {score}\n")
        f.write(f"Genome: {a1}\n")
        f.write(f"Gene:   {a2}\n")
        f.write("-------------------------\n")

res.sort(key=lambda x: x[0])
top10 = res[:10]
with open("top10Alignments.txt", "w") as f:
    for score, pos, a1, a2 in top10:
        f.write(f"Seed {pos}  score {score}\n")
        f.write(f"Genome: {a1}\n")
        f.write(f"Gene:   {a2}\n")
        f.write("-------------------------\n")