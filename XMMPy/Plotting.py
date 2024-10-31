import pandas as pd
import numpy as np
from XMMPy import hr_values

def rate_plotter(plt, file, band, error = "No", c = "blue", dashed = "No"):
    """Plot Count Rate Data

    Args:
        plt (MatPlotLib Pyplot): MatPlotLib Pyplot
        file (str): Absolute path to Count Rate file
        band (str): Name of energy band to be plotted
        error (str, optional): Whether to plot error bars. Defaults to "No".
        c (str, optional): Color of graph. Defaults to "blue".
        dashed (str, optional): Whether to make line dashed or not. Defaults to "No".
    """

    rates = pd.read_csv(file)
    xdata = (rates["Time"] - min(rates["Time"]))/1000
    ydata = rates[band]

    #Make dashed plot or not
    if dashed == "Yes":
        plt.step(xdata, ydata, where = "post", color = c, linestyle = "dotted", label = band)
    else:
        #Plot error bars
        if error == "Yes":
            errors = rates["Error"]
            plt.errorbar(xdata, ydata, yerr = errors, fmt = "o", color = "black", ecolor = "black", elinewidth = 0.75, capsize = 2, capthick = 0.75)
            plt.scatter(xdata, ydata, color = "black")
            running_average_plotter(plt, file, xdata, ydata)
        else:
            plt.step(xdata, ydata, where = "post", color = c, label = band)

    #Label axes
    plt.xlabel("Time (ks)", fontsize = 9)
    plt.ylabel("Count Rate (Cts/s)", fontsize = 9)

def counts_plotter(plt, file, band, c = "blue", dashed = "No"):
    """Plot Counts Data

    Args:
        plt (MatPlotLib Pyplot): MatPlotLib Pyplot
        file (str): Absolute path to Counts file
        band (str): Name of energy band to be plotted
        c (str, optional): Color of graph. Defaults to "blue".
        dashed (str, optional): Whether to make line dashed or not. Defaults to "No".
    """

    counts = pd.read_csv(file)
    xdata = (counts["Time"] - min(counts["Time"]))/1000
    ydata = counts[band]

    #Make dashed plot or not
    if dashed == "Yes":
        plt.step(xdata, ydata, where = "post", color = c, linestyle = "dotted", label = band)
    elif c == "black" and dashed == "No":
        plt.step(xdata, ydata, where = "post", color = c, label = band)
    else:
        plt.step(xdata, ydata, where = "post", color = c, label = band)

    #Label axes
    plt.xlabel("Time (ks)", fontsize = 9)
    plt.ylabel("Counts (Cts)", fontsize = 9)

def cumulative_counts_plotter(plt, file):
    """Plot Cumulative Counts

    Args:
        plt (MatPlotLib Pyplot): MatPlotLib Pyplot
        file (str): Absolute path to Counts file
    """

    counts = pd.read_csv(file)   
    xdata = (counts["Time"] - min(counts["Time"]))/1000
    ydata = counts["Broadband"]
    cumulative_counts = []
    running_total = 0

    #Accumulating values
    for count in ydata:
        running_total += count
        cumulative_counts.append(running_total)

    #Calculate values
    cumulative_counts = np.array(cumulative_counts)
    cumulative_counts = np.where(cumulative_counts == 0, np.nan, cumulative_counts)

    #Plot data and label axes
    plt.plot(xdata, cumulative_counts, color = "m")
    plt.xlabel("Time (ks)", fontsize = 9)
    plt.ylabel("Counts (Cts)", fontsize = 9)

def hr_plotter(plt, file):
    """Plot Hardness Ratios

    Args:
        plt (MatPlotLib Pyplot): MatPlotLib Pyplot
        file (str): Absolute path to Count Rate file
    """

    rates = pd.read_csv(file)
    xdata = (rates["Time"] - min(rates["Time"]))/1000
    for _, row in hr_values.iterrows():
        if row["Identifier"] == "m-s":
            ydata = ((rates["Medium"] - rates["Soft"]) / (rates["Medium"] + rates["Soft"]))
            plt.step(xdata, ydata, where = "post", color = row["Color"])
        elif row["Identifier"] == "h-s":
            ydata = ((rates["Hard"] - rates["Soft"]) / (rates["Hard"] + rates["Soft"]))
            plt.step(xdata, ydata, where = "post", color = row["Color"])
        elif row["Identifier"] == "s-m-h":
            ydata = (((rates["Soft"] - (rates["Medium"] + rates["Hard"])) / (rates["Soft"] + (rates["Medium"] + rates["Hard"]))))
            plt.step(xdata, ydata, where = "post", color = row["Color"])
        else:
            ydata = ((((rates["Hard"] + rates["Medium"]) - (rates["Soft"] + rates["Ultrasoft"])) / ((rates["Hard"] + rates["Medium"]) + (rates["Soft"] + rates["Ultrasoft"]))))
            plt.step(xdata, ydata, where = "post", color = row["Color"])

    #Label axes
    plt.xlabel("Time (ks)", fontsize = 9)
    plt.ylabel("Hardness Ratio", fontsize = 9)

def running_average_plotter(plt, file, xdata, ydata, window = 2):
    """Plot running average

    Args:
        plt (MatPlotLib Pyplot): MatPlotLib Pyplot
        file (str): Absolute path to Count Rate or Counts file
        xdata (Pandas Series): Series of "Time"
        ydata (Pandas Series): Series of "Rate" or "Counts"
        window (int, optional): Window size of running average. Defaults to 2.
    """
    
    df = pd.read_csv(file)
    averages = []
    num_of_bins = len(df["Bin"])
    for bin in df["Bin"]:
        index = bin - 1
        if index < 0 or index >= num_of_bins:
            averages.append(ydata[index] if 0 <= index < len(ydata) else np.nan)
        else:
            start_index = max(0, index - window)
            end_index = min(num_of_bins, index + window + 1)
            averages.append(np.mean(ydata[start_index:end_index]))

    plt.plot(xdata, averages, color = "red")