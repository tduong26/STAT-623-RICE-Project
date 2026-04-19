import random
G = 4000
L_min, L_max = 750, 800
L_mean = 775
N = 154
trials = 1000

gaps = 0
for _ in range(trials):
    reads = []
    for _ in range(N):
        l = random.randint(L_min, L_max)
        s = random.randint(0, G - l)
        reads.append({'start': s, 'end': s + l})
    
    sorted_reads = sorted(reads, key=lambda x: x['start'])
    if not sorted_reads:
        pass
    else:
        c_end = sorted_reads[0]['end']
        contigs = 1
        for r in sorted_reads[1:]:
            if r['start'] <= c_end:
                c_end = max(c_end, r['end'])
            else:
                contigs += 1
                c_end = r['end']
        if contigs > 1:
            gaps += contigs - 1

print("Avg contigs:", 1 + gaps/trials)
