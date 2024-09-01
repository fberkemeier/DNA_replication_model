import numpy as np
import os

chr_lengths = [249251, 243200, 198023, 191155, 180916, 171116, 159139, 146365, 141214, 135535, 135007, 133852, 115170, 107350, 102532, 90355, 81196, 78078, 59129, 63026, 48130, 51305]

### MAIN PROCESS ###

def process_bcs_output(cell_line, chr_number, chrpos_min, chrpos_max, fork_speed, resolution, scale_factor, sim_number, compute_replication_time, compute_fork_directionality, compute_origin_positions):
    # Define the file path
    file_path = f'data/bcs_output/bcs_output_{cell_line}_chr[{chr_number}]_{chrpos_min}-{chrpos_max}.simulation.bcs'

    # Initialize arrays to store replication time and fork directionality
    DNA_replicationtime = [0.0 for _ in range(0, chrpos_max - chrpos_min)] if compute_replication_time else None
    DNA_forkdirectionality = [0.0 for _ in range(0, chrpos_max - chrpos_min)] if compute_fork_directionality else None
    DNA_originpositions = [] if compute_origin_positions else None  # List to store origin positions per simulation
    current_origins = []
    sim_iteration = 0

    with open(file_path) as f:
        for line in f:
            if sim_iteration == sim_number + 1:
                break
            if line[0] == '>':
                alreadyDone = []
                if compute_origin_positions and current_origins:  # If we have collected origins for the current simulation
                    DNA_originpositions.append(current_origins)
                current_origins = []
                print(sim_iteration, end="\r")
                sim_iteration += 1
                continue
            splitLine = line.split('\t')
            if compute_origin_positions and splitLine[2] == "ORI":
                origin_pos = int(splitLine[4])
                current_origins.append(origin_pos)
            if splitLine[2] == "FL":
                pos = int(splitLine[4]) - 1
                time = float(splitLine[0])
                if pos not in alreadyDone:
                    if compute_replication_time:
                        DNA_replicationtime[pos] += time
                    if compute_fork_directionality:
                        DNA_forkdirectionality[pos] -= 1  # Track left-moving forks
                    alreadyDone.append(pos)
            if splitLine[2] == "FR":
                pos = int(splitLine[4]) - 1
                time = float(splitLine[0])
                if pos not in alreadyDone:
                    if compute_replication_time:
                        DNA_replicationtime[pos] += time
                    if compute_fork_directionality:
                        DNA_forkdirectionality[pos] += 1  # Track right-moving forks
                    alreadyDone.append(pos)

    # Don't forget to add the origins of the last simulation
    if compute_origin_positions and current_origins:
        DNA_originpositions.append(current_origins)

    # Average the results over the number of simulations
    if compute_replication_time:
        for i in range(len(DNA_replicationtime)):
            DNA_replicationtime[i] = float(DNA_replicationtime[i]) / float(sim_number)

    if compute_fork_directionality:
        for i in range(len(DNA_forkdirectionality)):
            DNA_forkdirectionality[i] = float(DNA_forkdirectionality[i]) / float(sim_number)

    # Define file paths for saving the results
    base_path = 'data'
    replication_time_path = os.path.join(base_path, 'whole-genome_timing_bcs', f'time_bcs_{cell_line}_chr[{chr_number}]_{chrpos_min}-{chrpos_max}.txt')
    fork_directionality_path = os.path.join(base_path, 'whole-genome_fork_directionality', f'fork_directionality_{cell_line}_chr[{chr_number}]_{chrpos_min}-{chrpos_max}.txt')
    origin_positions_path = os.path.join(base_path, 'whole-genome_origins', f'origin_positions_{cell_line}_chr[{chr_number}]_{chrpos_min}-{chrpos_max}.txt')

    # Create directories if they do not exist
    os.makedirs(os.path.dirname(replication_time_path), exist_ok=True)
    os.makedirs(os.path.dirname(fork_directionality_path), exist_ok=True)
    os.makedirs(os.path.dirname(origin_positions_path), exist_ok=True)

    # Save the results to text files
    if compute_replication_time:
        np.savetxt(replication_time_path, DNA_replicationtime, fmt='%.6f')

    if compute_fork_directionality:
        np.savetxt(fork_directionality_path, DNA_forkdirectionality, fmt='%.6f')

    if compute_origin_positions:
        with open(origin_positions_path, 'w') as f:
            for origins in DNA_originpositions:
                f.write(' '.join(map(str, origins)) + '\n')



### PROCESS INTERVALS ###

def process_intervals(cell_lines, chr_numbers, fork_speed=1.4, resolution=1000, scale_factor=6, sim_number=5, compute_replication_time=True, compute_fork_directionality=True, compute_origin_positions=True, interval=None):
    for cell_line in cell_lines:
        for chr_number in chr_numbers:
            chr_length = chr_lengths[chr_number - 1]  # Get the length of the chromosome
            intervals = [(interval[0], interval[1])] if interval else [(start, min(start + 10000, chr_length)) for start in range(0, chr_length, 10000)]
            for start, end in intervals:
                process_bcs_output(
                    cell_line=cell_line,
                    chr_number=chr_number,
                    chrpos_min=start,
                    chrpos_max=end,
                    fork_speed=fork_speed,
                    resolution=resolution,
                    scale_factor=scale_factor,
                    sim_number=sim_number,
                    compute_replication_time=compute_replication_time,
                    compute_fork_directionality=compute_fork_directionality,
                    compute_origin_positions=compute_origin_positions
                )
                
### INTERORIGIN DISTANCE ###

def compute_interorigin_distances(cell_line, chr_number, chrpos_min, chrpos_max):
    base_path = 'data'
    
    # Define file path for loading the origins
    origins_path = os.path.join(base_path, 'whole-genome_origins', f'origin_positions_{cell_line}_chr[{chr_number}]_{chrpos_min}-{chrpos_max}.txt')
    
    # Load origins data from text file
    if os.path.exists(origins_path):
        with open(origins_path, 'r') as f:
            origins_data = [list(map(int, line.strip().strip('[]').split())) for line in f]
    else:
        raise FileNotFoundError(f"Origins data not found at {origins_path}")
    
    # Compute interorigin distances for each simulation
    interorigin_distances = []
    for origins in origins_data:
        origins_sorted = sorted(origins)
        distances = np.diff(origins_sorted)
        interorigin_distances.append(distances)
    
    # Define file path for saving the interorigin distances
    iod_path = os.path.join(base_path, 'whole-genome_interorigin_distances', f'iod_{cell_line}_chr[{chr_number}]_{chrpos_min}-{chrpos_max}.txt')
    
    # Save interorigin distances to a text file
    with open(iod_path, 'w') as f:
        for distances in interorigin_distances:
            f.write(f"{list(distances)}\n")

def compute_interorigin_intervals(cell_lines, chr_numbers, interval=None):
    for cell_line in cell_lines:
        for chr_number in chr_numbers:
            chr_length = chr_lengths[chr_number - 1]  # Get the length of the chromosome
            intervals = [(interval[0], interval[1])] if interval else [(start, min(start + 10000, chr_length)) for start in range(0, chr_length, 10000)]
            for start, end in intervals:
                compute_interorigin_distances(
                    cell_line=cell_line,
                    chr_number=chr_number,
                    chrpos_min=start,
                    chrpos_max=end
                )

def average_iod_data(cell_lines, chr_numbers, factor_min=5, show_per_cell_line=False):
    all_iod_data = []
    iod_data_per_cell_line = []

    for cell_line in cell_lines:
        cell_line_iod_data = []
        for chr_number in chr_numbers:
            chr_length = chr_lengths[chr_number - 1]
            for start in range(0, chr_length, 10000):
                end = min(start + 10000, chr_length)
                iod_data = load_function_metrics(cell_line, chr_number, "iod", start, end, factor_min=factor_min)
                cell_line_iod_data.extend(iod_data)
        iod_data_per_cell_line.append(cell_line_iod_data)
        all_iod_data.extend(cell_line_iod_data)

    if show_per_cell_line:
        return iod_data_per_cell_line
    else:
        return [all_iod_data]
        
### EXAMPLE USAGE ###

# Example usage
cell_lines = ["HCT"]
chr_numbers = [1]
sim_number=500

process_intervals(cell_lines, chr_numbers, sim_number=sim_number)