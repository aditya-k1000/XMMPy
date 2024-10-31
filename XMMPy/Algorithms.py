import pandas as pd

def eclipse_detect(binsize, counts_file):
    """Checks for eclipses in counts file

    Args:
        binsize (int): Binsize (s)
        counts_file (str): Absolute path to counts file

    Returns:
        bool: Boolean value based on presence of eclipse or not
    """
    
    df = pd.read_csv(counts_file)

    if len(df) < 5:
        return False

    broadband_std = df["Broadband"].std()
    threshold = 0.5 * broadband_std
    dip_threshold = 10 * (binsize / 500)
    consecutive_count = 0

    for i in range(len(df) - 1):
        if abs(df["Broadband"][i] - df["Broadband"][i + 1]) >= dip_threshold:
            consecutive_count = 1
            for j in range(i + 1, len(df) - 1):
                if abs(df["Broadband"][j] - df["Broadband"][j + 1]) <= threshold:
                    consecutive_count += 1
                    if consecutive_count >= 3:
                        return True
                else:
                    break
    return False