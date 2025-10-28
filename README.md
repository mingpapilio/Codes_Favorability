# Codes_Favorability

This repository has the codes and data for the manuscript "__Warmer environments harbor greater thermal trait diversity__".  There are two sections in this repository: simulations and empirical analysis. 

## Description of the Data and File Structure

The repository is structured into four main directories, each corresponding to a major analysis section of the manuscript:

### 1. **`Simulations`**
This folder contains the (1) source codes for the eco-evolutionary individual-based model, and (2) the simulated data used in making main and extended figures.

- `codes`: 

- `dSFMT-src-2.3.3`: This is a dependency for the standard random number generator (Fast Mersenne Twister), equivalent to std::mt19937 in C++. Source code is available here `https://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/SFMT/`.

- `data`: The simulation results are here for reproduciability

### 2. **`Empirical analysis`**
This folder contains our empirical data and code that are able to reproduce the results of Figure 4 and 5.

- `FAV_hv_data.csv`: data file that contains moth species level traits and corresponding ambient temperatures.
    - Coulmn names:
        - ``Elev_up``: Upper distribution limit
        - ``Elev_low``: Lower distribution limit
        - ``Elev_mid``: Mid-point of upper and lower distribution limits
        - ``RS``: Range size (Elev_up - Elev_low)
        - ``W_length``: Forewing length
        - ``B_length``: Body length
        - ``Weight``: Dry weight
        - ``CTmax``, ``CTmin``: Critical thermal maximum / minimum
        - ``TTrange``: Thermal tolerance range (CTmax - CTmin)
        - ``STmax``: Average max temperature of the warmest month
        - ``STmin``: Average min temperature of the coldest month
        - ``STR``: Seasonal temperature range (STmax - STmin)
        - ``STmean``: Average mean annual temperature
        - ``Location``: Location where the species was collected
        - ``DTmax``: Diurnal maximum temperature
        - ``DTmin``: Diurnal minimum temperature
        - ``DTR``: Diurnal temperature range (DTmax - DTmin)
        - ``DTmean``: Diurnal mean temperature
        - ``Group``: Group which species correspond to when applying assemblage-level analyses. For the meaning of the codes, please refer to Figure 4(a) in the article.

- `/Fig4/`: This folder contains all codes needed to produce Figure 4b-j.
    - R scripts:
        - ``Fig4_hypervolume.R``: For hypervolume calculation and plotting in Figure 4b, 4e, 4h.
        - ``Fig4_TTrange_histogram.R``: For plotting histograms in Figure 4c, 4f, 4i.
        - ``Fig4_TTrange_MeanSD.R``: For calculating mean and SD of TTrange of each assemblage and plotting in Figure 4d, 4g, 4j.
    - CSV files:
        - ``hv_result.csv``: The hypervolume values used for analysis in this article, which should be reproducible by using 'Fig4_hypervolume.R'.
            - Column names:
                - ``Location``: Location of the assemblage.
                - ``Group``: Group code which species correspond to when applying assemblage-level analyses. For the meaning of the codes, please refer to Figure 4(a) in the article.
                - ``Thermal_hv``: Average hypervolume of thermal tolerance traits (CTmax, CTmin) of assemblage after 100 iterations.
                - ``Tmean_all``: Averange mean ambient temperature of assemblage.
                - ``DTR_all``: Averange diurnal temperature range of assemblage.
                - ``STR_all``: Averange seasonal temperature range of assemblage.
                - ``group_num``: Contains only the numeric portion of the Group identifier (e.g., "M1" becomes 1, "T3" becomes 3). 
                - ``N_all``: Number of species of assemblage.


If you have any questions, please contact *ming.liu.ac[@]gmail.com* for simulation-related enquiries and *tzumanhung9527[@]gmail.com* for empirical-related enquiries.
