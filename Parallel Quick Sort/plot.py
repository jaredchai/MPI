import pandas as pd
import matplotlib.pyplot as plt

# CSV data as a string
csv_data = "results.csv"

# Read the CSV data using pandas
data = pd.read_csv(csv_data, header=0, index_col=0)

# Set the index (n) and columns (p) names
data.index.name = "n"
data.columns.name = "p"

# Plot the data
plt.figure(figsize=(10, 6))
for col in data.columns:
    plt.plot(data.index, data[col], marker='o', label=f'p = {col}')
plt.xlabel("Input size (n)")
plt.ylabel("Runtime (s)")
plt.title("Runtime vs Input Size for Different Number of Processors")
plt.xscale('log')
plt.legend()
plt.grid()
plt.show()


# Transpose the DataFrame
data = data.T

# Set the index (p) and columns (n) names
data.index.name = "p"
data.columns.name = "n"

# Plot the data
plt.figure(figsize=(10, 6))
for col in data.columns:
    plt.plot(data.index, data[col], marker='o', label=f'n = {col}')
plt.xlabel("Number of Processors (p)")
plt.ylabel("Runtime (s)")
plt.title("Runtime vs Number of Processors for Different Input Sizes")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.grid()
plt.tight_layout()
plt.show()
