import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from classes.zero_temp import Chain


for j in range(10):
    system1 = Chain(100000 ,100)

    for i in tqdm(range(10000)):
        system1.elimination_transformation()

print(len(system1.bonds))
