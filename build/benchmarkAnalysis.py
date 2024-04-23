import statistics

# Open the file
file0 = 'slurm-516331.out' # Sphere Wave Fixed
file1 = 'slurm-516332.out' # Sphere Wave Adaptive
file2 = 'slurm-516333.out' # Moving Diffusion Fixed
file3 = 'slurm-516334.out' # Moving Diffusion Adaptive
file4 = 'slurm-516335.out' # Curved Beam Fixed
file5 = 'slurm-516336.out' # Curved Beam Adaptive

with open(file5, 'r') as file:
    moments = []
    stream = []
    collide = []
    computation = []

    for line in file:
        if 'ComputeMomentsIF' in line:
            value = line.split('():')[1].strip()
            moments.append(float(value[:-1]))
        if 'Stream' in line:
            value = line.split('():')[1].strip()
            stream.append(float(value[:-1]))
        if 'Collide' in line:
            value = line.split('():')[1].strip()
            collide.append(float(value[:-1]))
        if 'Computation Time' in line:
            value = line.split(':')[1].strip()
            computation.append(float(value[:-1]))

print("Moments Average:", statistics.mean(moments))
print("Moments Std Dev:", statistics.stdev(moments))
print("Stream  Average:", statistics.mean(stream))
print("Stream  Std Dev:", statistics.stdev(stream))
print("Collide Average:", statistics.mean(collide))
print("Collide Std Dev:", statistics.stdev(collide))
print("Total   Average:", statistics.mean(computation))
print("Total   Std Dev:", statistics.stdev(computation))