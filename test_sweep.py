"""
test_sweep.py
=============
Standalone read-length sweep test.
Demonstrates how perfect (snp_rate=0) repeats affect assembly fragmentation.
No pysam dependency.
"""
import shortgun_simulation as ss

G_sweep = 15000
repeat_sizes_sweep = {"LINE": 800, "SINE": 300, "LTR": 400}
weights_sweep = [
    ("LINE",   0.25),
    ("SINE",   0.15),
    ("LTR",    0.16),
    ("UNIQUE", 0.44),
]

# snp_rate=0.0 → repeats are perfectly identical, worst-case for assembly
genome_sweep = ss.generate_genome_with_repeats(
    G_sweep,
    snp_rate=0.00,
    repeat_sizes=repeat_sizes_sweep,
    weights=weights_sweep,
)

ss.run_read_length_sweep(
    genome_sweep,
    target_coverage=30,
    output_png="read_length_sweep_s2.png",
    label="Perfect Repeats",
    num_trials=3,
)
