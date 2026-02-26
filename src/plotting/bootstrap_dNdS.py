from plot_dNdS import get_summary_paths,get_species_list,violinplot_pair,violinplot_pair_single
from plot_dS import read_dNdS_dS_summary_file
import scipy.stats as sts
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter,FuncFormatter
import analyze_site_classes as site_classes



def read_filtered_dNdS_summary(summary_path, excl_list=[], max_dS=2, max_dNdS = 2):
    """
    read the dNdS and dS values, filter for min_dS, and return only dNdS values that meet the criteria
    """
    summary_dict = read_dNdS_dS_summary_file(summary_path=summary_path, exclude_list=excl_list, only_dS=False)
    dNdS_dict = {pair : [] for pair in summary_dict.keys()}
    for pair, lists_dict in summary_dict.items():
        dNdS_unfiltered = lists_dict["dNdS"]
        dS_list = lists_dict["dS"]
        assert len(dNdS_unfiltered) == len(dS_list)

        dNdS_filtered = [dNdS for i, dNdS in enumerate(dNdS_unfiltered) if dS_list[i]<max_dS and dNdS<max_dNdS]
        print(f"{pair} : {len(dNdS_unfiltered)} dNdS values, {len(dS_list)} dS values; \t {len(dNdS_filtered)} have dS < {max_dS}")
        dNdS_dict[pair] = dNdS_filtered

    return dNdS_dict


def permute_dNdS(dNdS_A, dNdS_X):
    """
    return resampled A and X
    """
    n_A = len(dNdS_A)
    dNdS_all = dNdS_A + dNdS_X
    permut = np.random.permutation(dNdS_all) ## permutation: sampling without replacement
    new_A = permut[0:n_A]
    new_X = permut[n_A:]
    return new_A, new_X


def permutate_dNdS(dNdS_A, dNdS_X, num_permut = 1000, mean=False):
    """
    permutate n times and calculate the differences of median between all pairs
    """
    medians_diff_list = [np.NaN] * num_permut
    print(f"... running {num_permut} permutations ...")
    if mean:
        for i in range(num_permut):
            new_A, new_X = permute_dNdS(dNdS_A=dNdS_A, dNdS_X=dNdS_X)
            m_A = np.nanmedian(new_A)
            m_X = np.nanmedian(new_X)
            medians_diff_list[i] = m_A - m_X
    else: # median
        for i in range(num_permut):
            new_A, new_X = permute_dNdS(dNdS_A=dNdS_A, dNdS_X=dNdS_X)
            m_A = np.nanmean(new_A)
            m_X = np.nanmean(new_X)
            medians_diff_list[i] = m_A - m_X
    
    return medians_diff_list



def calculate_list_CI(values_list:list, cl = 0.95, verbose = False):
    """
    calculate 95% confidence interval of a list of float values
    """
    mean_coeff = np.mean(values_list)
    std_coeff = np.std(values_list)
            
    # Sample statistics
    lower, upper = sts.norm.interval(cl, loc = mean_coeff, scale = std_coeff) 
    norm_coeffs = [mean_coeff,std_coeff, lower, upper]
    if verbose:
        print(f"\t\tmean correlation coefficient: {mean_coeff:.3f}, standard deviation {std_coeff:.3f}, 95% confidence interval: [{lower:.3f}, {upper:.3f}]")
    # ci = sts.t.interval(cl, df=len(values_list)-1, loc=np.mean(values_list), scale=np.std(values_list, ddof=1) / np.sqrt(len(values_list)))
    return(norm_coeffs)



def plot_dNdS_permutations(boot_diff:dict, measure_diff:dict, A_dict:dict, X_dict:dict, filename = "dNdS_permutations.png", dark_mode=False, hist_label = r"$\text{dNdS}_A - \text{dNdS}_X$", violin_label="dNdS", violin_ymax = 0, transparent=True, binary = False):
    """
    plot a grid of histogram distribution plots for all pairwise comparisons
    """

    if dark_mode:
        plt.style.use('dark_background')

    species_list = get_species_list(boot_diff)
    print(f"... plotting {len(species_list)} species")

    species_count = len(species_list)
    species_index = {species : i for i, species in enumerate(species_list)}

    cols = species_count
    rows = cols
    if rows>2:
        fig, axes = plt.subplots(rows, cols, figsize=(27, 25)) # for more than three rows
    else:
        fig, axes = plt.subplots(rows, cols, figsize=(15, 10)) # for less than three rows
    
    fs = 30

    plt.rcParams['text.usetex'] = True

    colors_dict = {
        "grey" : "#3A4040", 
        "bars" : "#7BB7AE", 
        "line" : "#FE4894",
        "pdf" : "#3D7068", 
        "A" : "#F2933A", # uniform_filtered orange
        "X" : "#b82946", # native red
    }

    diagonals_done = []

    # sort so that the order stays the same every time. otherwise it changes
    pairs_sorted = sorted(list(boot_diff.keys()))
    for pair in pairs_sorted:
        ### get pair indices for species pair
        try:
            gen1, spec1, gen2, spec2 =pair.split("_")
        except:
            raise RuntimeError(f"{pair} could not be parsed")
        species1 = f"{gen1}_{spec1}"
        row = species_index[species1]
        species2 = f"{gen2}_{spec2}"
        col = species_index[species2]

        ## plot center row/col species names
        # put this here before otherwise the last row/col label never gets reached
        if row == len(species_list)-1:
            species1_lab = species1.replace("_", ". ")
            axes[row,row].text(0.1,0.4,f"{species1_lab}", fontsize = fs*1.4)
        if col == len(species_list)-1:
            species2_lab = species2.replace("_", ". ")
            axes[col,col].text(0.1,0.4,f"{species2_lab}", fontsize = fs*1.4)

        # only do top right matrix
        if row>col:
            col_temp = col
            col = row
            row = col_temp
            species2_temp = species2
            species2 = species1
            species1 = species2_temp
        species1_lab = species1.replace("_", ". ")
        species2_lab = species2.replace("_", ". ")

        ## plot species name on diagonals
        if row not in diagonals_done:
            axes[row,row+1].set_xlabel(hist_label, fontsize = fs)
            # axes[row,row].text(0.8,0.2,f"{species1_lab}", rotation = 90, fontsize = fs*1.3)
            axes[row,row].text(0.1,0.4,f"{species1_lab}", fontsize = fs*1.4)
            diagonals_done.append(row)
        
        ### plot bootstraps

        data_boot = np.array(boot_diff[pair], dtype=float)
    
        if len(data_boot)==0:
            axes[row,col].axis('off')
            axes[row,col].set_title(f'{species1_lab}\n{species2_lab}', fontsize = fs*0.85)
            continue

        n_boot = len(data_boot)
        mean_boot = np.nanmedian(data_boot)

        nbins = 20
        axes[row,col].hist(data_boot, bins=nbins, weights=np.ones(len(data_boot)) / len(data_boot), histtype="bar", color = [colors_dict["bars"]])#, density = True)
        axes[row,col].yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 99 and x<1 else f'{int(x*100)}%'))
        axes[row,col].axvline(measure_diff[pair], ymin = 0, ymax = len(data_boot)/nbins, color = colors_dict["line"], linewidth = 3, label = f"median(A) - median(X)")
        # axes[row,col].legend(fontsize = fs*0.85)
        # axes[row,col].axvline(0, ymin = 0, ymax = len(data_boot)/nbins, color = colors_dict["grey"], linewidth = 1, linestyle = "--")

        ## plot normal pdf based on sampling
        lwd = 3
        pdf_ax = axes[row,col].twinx()
        if measure_diff[pair] > max(data_boot):
            pdf_x = np.arange(min(data_boot),measure_diff[pair], 0.00001)
        else:
            pdf_x = np.arange(min(data_boot),max(data_boot), 0.00001)
        mean_cor,std_cor,lower_CI,upper_CI = calculate_list_CI(data_boot)
        pdf_y = sts.norm.pdf(pdf_x, mean_cor, std_cor)
        pdf_ax.plot(pdf_x,pdf_y, linewidth=lwd, color=colors_dict["pdf"], label = r"mean $\text{dNdS}_A - \text{dNdS}_X$"+f": {mean_cor:.3f}")#, linestyle="--")
        pdf_ax.axvline(x=lower_CI, linewidth=lwd, color=colors_dict["pdf"], linestyle=":", label = f"95% standard error \nof the mean [{lower_CI:.3f}, {upper_CI:.3f}]")
        pdf_ax.axvline(x=upper_CI, linewidth=lwd, color=colors_dict["pdf"], linestyle=":")
        axes[row,col].axvspan(lower_CI, upper_CI, color=colors_dict["pdf"], alpha=0.25)
        pdf_ax.axis('off')

        # axes[row,col].set_title(f'{species2_lab}', fontsize = fs)
        # axes[row,col].set_title(f'{species2}\n{species1}', fontsize = fs*0.85)
        axes[row,col].tick_params(axis='y', labelsize=fs)
        axes[row,col].tick_params(axis='x', labelsize=fs) 

        axes[row,row].axis('off')
        axes[col,col].axis('off')

        
        ### plot violins or binary bar charts

        ## exclude all the NaNs because violinplot can't handle them
        data_A_nan = np.array(A_dict[pair], dtype=float)
        data_X_nan = np.array(X_dict[pair], dtype=float)
        data_A = [dNdS_A for dNdS_A in data_A_nan if not np.isnan(dNdS_A) ]
        data_X = [dNdS_X for dNdS_X in data_X_nan if not np.isnan(dNdS_X) ]
    
        if len(data_A)==0 or len(data_X)==0:
            # axes[col,row].axis('off')
            # axes[col,row].set_title(f'{species1}\n{species2}', fontsize = fs*0.85)
            continue

        n_A = len(data_A)
        n_X = len(data_X)
        mean_A = np.nanmedian(data_A)
        mean_X = np.nanmedian(data_X)

        data_AX = [data_A, data_X]

        if binary:
            site_classes.binary_barplot_pair(data_A=data_A, data_X=data_X, row=col, col=row, n_A=n_A, n_X=n_X, axes=axes, colors_dict=colors_dict, fs = fs, ylab ="pos. sel. genes")
        else:
            violinplot_pair(data_A_X=data_AX, row=col, col=row, n_A=n_A, n_X=n_X, mean_A=mean_A, mean_X=mean_X, axes = axes, colors_dict=colors_dict, fs = fs, xlab = violin_label, ymax = violin_ymax)
        # axes[col,row].set_title(f'{species2}\n{species1}', fontsize = fs*0.85)
        # axes[col,row].set_title(f'{species2}', fontsize = fs)

        print(f"{row}, {col} : {species1_lab} vs. {species2_lab} --> mean(dNdS_A)-mean(dNdS_X) bootstrap: {mean_boot:.5f}, measured: {measure_diff[pair]:.3f}")

    fig.suptitle(f"Bruchini: A and X dNdS violin plot and permutation test\n ", fontsize = fs*1.5)
    # Adjust layout to prevent overlap  (left, bottom, right, top)
    plt.tight_layout(rect=[0.01, 0, 1, 1])

    if dark_mode:
        filename = filename.replace(".png", "_darkmode.png")

    # transparent background
    plt.savefig(filename, dpi = 300, transparent = True)
    # non-transparent background
    filename_tr = filename.replace(".png", "_white_bg.png")
    plt.savefig(filename_tr, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as: {filename} and {filename_tr}")


def plot_dNdS_permutations_one_pair(boot_diff:dict, measure_diff:dict, A_dict:dict, X_dict:dict, filename = "dNdS_permutations.png", dark_mode=False, hist_label = r"$\text{dNdS}_A - \text{dNdS}_X$", violin_label="dNdS", violin_ymax = 0, transparent=True, binary = False):
    """
    plot a grid of histogram distribution plots for all pairwise comparisons
    """

    if dark_mode:
        plt.style.use('dark_background')

    species_list = get_species_list(boot_diff)
    print(f"... plotting {len(species_list)} species")

    species_count = len(species_list)
    assert species_count==2
    species_index = {species : i for i, species in enumerate(species_list)}

    cols = species_count
    rows = 1
    fig, axes = plt.subplots(rows, cols, figsize=(16, 8)) # for less than three rows
    fs = 30
    plt.rcParams['text.usetex'] = True

    colors_dict = {
        "grey" : "#3A4040", 
        "bars" : "#7BB7AE", 
        "line" : "#FE4894",
        "pdf" : "#3D7068", 
        "A" : "#F2933A", # uniform_filtered orange
        "X" : "#b82946", # native red
    }

    pairs_sorted = sorted(list(boot_diff.keys()))
    assert len(pairs_sorted) ==1
    for pair in pairs_sorted:
        ### get pair indices for species pair
        try:
            gen1, spec1, gen2, spec2 =pair.split("_")
        except:
            raise RuntimeError(f"{pair} could not be parsed")
        species1 = f"{gen1}_{spec1}"
        species2 = f"{gen2}_{spec2}"

        species1_lab = species1.replace("_", ". ")
        species2_lab = species2.replace("_", ". ")
        fig.suptitle(f"{species1_lab} vs. {species2_lab}: \nA and X dNdS violin plot and permutation test", fontsize = fs*1.5)
        
        ### plot bootstraps

        data_boot = np.array(boot_diff[pair], dtype=float)

        n_boot = len(data_boot)
        mean_boot = np.nanmedian(data_boot)

        nbins = 20
        axes[1].hist(data_boot, bins=nbins, weights=np.ones(len(data_boot)) / len(data_boot), histtype="bar", color = [colors_dict["bars"]])#, density = True)
        axes[1].yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 99 and x<1 else f'{int(x*100)}%'))
        axes[1].axvline(measure_diff[pair], ymin = 0, ymax = len(data_boot)/nbins, color = colors_dict["line"], linewidth = 3, label = f"median(A) - median(X)")
        # axes[1].legend(fontsize = fs*0.85)
        # axes[row,col].axvline(0, ymin = 0, ymax = len(data_boot)/nbins, color = colors_dict["grey"], linewidth = 1, linestyle = "--")

        ## plot normal pdf based on sampling
        lwd = 3
        pdf_ax = axes[1].twinx()
        if measure_diff[pair] > max(data_boot):
            pdf_x = np.arange(min(data_boot),measure_diff[pair], 0.00001)
        else:
            pdf_x = np.arange(min(data_boot),max(data_boot), 0.00001)
        mean_cor,std_cor,lower_CI,upper_CI = calculate_list_CI(data_boot)
        pdf_y = sts.norm.pdf(pdf_x, mean_cor, std_cor)
        pdf_ax.plot(pdf_x,pdf_y, linewidth=lwd, color=colors_dict["pdf"], label = r"mean $\text{dNdS}_A - \text{dNdS}_X$"+f": {mean_cor:.3f}")#, linestyle="--")
        pdf_ax.axvline(x=lower_CI, linewidth=lwd, color=colors_dict["pdf"], linestyle=":", label = f"95% standard error \nof the mean [{lower_CI:.3f}, {upper_CI:.3f}]")
        pdf_ax.axvline(x=upper_CI, linewidth=lwd, color=colors_dict["pdf"], linestyle=":")
        axes[1].axvspan(lower_CI, upper_CI, color=colors_dict["pdf"], alpha=0.25)
        axes[1].set_xlabel(hist_label, fontsize = fs)
        pdf_ax.axis('off')

        # axes[row,col].set_title(f'{species2_lab}', fontsize = fs)
        # axes[row,col].set_title(f'{species2}\n{species1}', fontsize = fs*0.85)
        axes[1].tick_params(axis='y', labelsize=fs)
        axes[1].tick_params(axis='x', labelsize=fs) 
        
        ### plot violins or binary bar charts

        ## exclude all the NaNs because violinplot can't handle them
        data_A_nan = np.array(A_dict[pair], dtype=float)
        data_X_nan = np.array(X_dict[pair], dtype=float)
        data_A = [dNdS_A for dNdS_A in data_A_nan if not np.isnan(dNdS_A) ]
        data_X = [dNdS_X for dNdS_X in data_X_nan if not np.isnan(dNdS_X) ]
    
        if len(data_A)==0 or len(data_X)==0:
            # axes[col,row].axis('off')
            # axes[col,row].set_title(f'{species1}\n{species2}', fontsize = fs*0.85)
            continue

        n_A = len(data_A)
        n_X = len(data_X)
        mean_A = np.nanmedian(data_A)
        mean_X = np.nanmedian(data_X)

        data_AX = [data_A, data_X]

        if binary:
            site_classes.binary_barplot_single_pair(data_A=data_A, data_X=data_X, row=0, n_A=n_A, n_X=n_X, axes=axes, colors_dict=colors_dict, fs = fs, ylab ="pos. sel. genes")
        else:
            violinplot_pair_single(data_A_X=data_AX, col=0, n_A=n_A, n_X=n_X, mean_A=mean_A, mean_X=mean_X, axes = axes, colors_dict=colors_dict, fs = fs, xlab = violin_label, ymax = violin_ymax)
        # axes[col,row].set_title(f'{species2}\n{species1}', fontsize = fs*0.85)
        # axes[col,row].set_title(f'{species2}', fontsize = fs)

        print(f"{species1_lab} vs. {species2_lab} --> mean(dNdS_A)-mean(dNdS_X) bootstrap: {mean_boot:.5f}, measured: {measure_diff[pair]:.3f}")

    # Adjust layout to prevent overlap (left, bottom, right, top)
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    if dark_mode:
        filename = filename.replace(".png", "_darkmode.png")

    # transparent background
    plt.savefig(filename, dpi = 300, transparent = True)
    # non-transparent background
    filename_tr = filename.replace(".png", "_white_bg.png")
    plt.savefig(filename_tr, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as: {filename} and {filename_tr}")


def plot_dNdS_pos_sel(A_dict_dNdS:dict, X_dict_dNdS:dict, pos_list_A:dict, pos_list_X:dict, filename = "dNdS_pos_sel.png", dark_mode=False, violin_label="dNdS", violin_ymax = 0):
    """
    plot a grid of histogram distribution plots for all pairwise comparisons
    """

    if dark_mode:
        plt.style.use('dark_background')

    species_list = get_species_list(A_dict_dNdS)
    print(f"... plotting {len(species_list)} species")

    species_count = len(species_list)
    species_index = {species : i for i, species in enumerate(species_list)}

    cols = species_count
    rows = cols
    if rows>2:
        fig, axes = plt.subplots(rows, cols, figsize=(27, 25)) # for more than three rows
    else:
        fig, axes = plt.subplots(rows, cols, figsize=(15, 10)) # for less than three rows
    
    fs = 30

    plt.rcParams['text.usetex'] = True

    colors_dict = {
        "grey" : "#3A4040", 
        "bars" : "#7BB7AE", 
        "line" : "#FE4894",
        "pdf" : "#3D7068", 
        "A" : "#F2933A", # uniform_filtered orange
        "X" : "#b82946", # native red
    }

    diagonals_done = []

    # sort so that the order stays the same every time. otherwise it changes
    pairs_sorted = sorted(list(A_dict_dNdS.keys()))
    for pair in pairs_sorted:

        ### get pair indices for species pair
        try:
            gen1, spec1, gen2, spec2 =pair.split("_")
        except:
            raise RuntimeError(f"{pair} could not be parsed")
        species1 = f"{gen1}_{spec1}"
        row = species_index[species1]
        species2 = f"{gen2}_{spec2}"
        col = species_index[species2]

        ## plot center row/col species names
        # put this here before otherwise the last row/col label never gets reached
        if row == len(species_list)-1:
            species1_lab = species1.replace("_", ". ")
            axes[row,row].text(0.1,0.4,f"{species1_lab}", fontsize = fs*1.4)
        if col == len(species_list)-1:
            species2_lab = species2.replace("_", ". ")
            axes[col,col].text(0.1,0.4,f"{species2_lab}", fontsize = fs*1.4)

        # only do top right matrix
        axes[row,row].axis('off')
        axes[col,col].axis('off')
        if row>col:
            col_temp = col
            col = row
            row = col_temp
            species2_temp = species2
            species2 = species1
            species1 = species2_temp
        species1_lab = species1.replace("_", ". ")
        species2_lab = species2.replace("_", ". ")


        ## plot species name on diagonals
        if row not in diagonals_done:
            # axes[row,row].text(0.8,0.2,f"{species1_lab}", rotation = 90, fontsize = fs*1.3)
            axes[row,row].text(0.1,0.4,f"{species1_lab}", fontsize = fs*1.4)
            diagonals_done.append(row)
    
        
        ### plot binary bar charts
        ## exclude all the NaNs because violinplot can't handle them
        data_A_nan = np.array(pos_list_A[pair], dtype=float)
        data_X_nan = np.array(pos_list_X[pair], dtype=float)
        data_A = [dNdS_A for dNdS_A in data_A_nan if not np.isnan(dNdS_A) ]
        data_X = [dNdS_X for dNdS_X in data_X_nan if not np.isnan(dNdS_X) ]
    
        if len(data_A)==0 or len(data_X)==0:
            # axes[col,row].axis('off')
            # axes[col,row].set_title(f'{species1}\n{species2}', fontsize = fs*0.85)
            continue

        n_A = len(data_A)
        n_X = len(data_X)
        mean_A = np.nanmedian(data_A)
        mean_X = np.nanmedian(data_X)

        data_AX = [data_A, data_X]

        print(f"{row}, {col} : {species1_lab} vs. {species2_lab} --> sig.pos.sel bar plot")
        
        site_classes.binary_barplot_pair(data_A=data_A, data_X=data_X, row=row, col=col, n_A=n_A, n_X=n_X, axes=axes, colors_dict=colors_dict, fs = fs, ylab ="pos. sel. genes")
        
        ### plot violins
        ## exclude all the NaNs because violinplot can't handle them
        try:
            data_A_nan = np.array(A_dict_dNdS[pair], dtype=float)
            data_X_nan = np.array(X_dict_dNdS[pair], dtype=float)
        except:
            raise RuntimeError(f"data does not have key: {pair}, available keys are:\n{A_dict_dNdS.keys()}")
        data_A = [dNdS_A for dNdS_A in data_A_nan if not np.isnan(dNdS_A) ]
        data_X = [dNdS_X for dNdS_X in data_X_nan if not np.isnan(dNdS_X) ]
    
        if len(data_A)==0 or len(data_X)==0:
            # axes[col,row].axis('off')
            # axes[col,row].set_title(f'{species1}\n{species2}', fontsize = fs*0.85)
            continue

        n_A = len(data_A)
        n_X = len(data_X)
        mean_A = np.nanmedian(data_A)
        mean_X = np.nanmedian(data_X)

        data_AX = [data_A, data_X]


        violinplot_pair(data_A_X=data_AX, row=col, col=row, n_A=n_A, n_X=n_X, mean_A=mean_A, mean_X=mean_X, axes = axes, colors_dict=colors_dict, fs = fs, xlab = violin_label, ymax = violin_ymax, text_x_offset=0.5)
        # axes[col,row].set_title(f'{species2}\n{species1}', fontsize = fs*0.85)
        # axes[col,row].set_title(f'{species2}', fontsize = fs)

        print(f"{col}, {row} : {species1_lab} vs. {species2_lab} --> dNdS violin plot\n")

    fig.suptitle(f"Bruchini: A and X dNdS violin plot and pos. selected genes\n ", fontsize = fs*1.5)
    # Adjust layout to prevent overlap  (left, bottom, right, top)
    plt.tight_layout(rect=[0.01, 0, 1, 1])

    if dark_mode:
        filename = filename.replace(".png", "_darkmode.png")

    # transparent background
    plt.savefig(filename, dpi = 300, transparent = True)
    # non-transparent background
    filename_tr = filename.replace(".png", "_white_bg.png")
    plt.savefig(filename_tr, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as: {filename} and {filename_tr}")



def plot_dNdS_pos_sel_one_pair(A_dict_dNdS:dict, X_dict_dNdS:dict, pos_list_A:dict, pos_list_X:dict, filename = "dNdS_pos_sel.png", dark_mode=False, violin_label="dNdS", violin_ymax = 0):
    """
    plot a grid of histogram distribution plots for all pairwise comparisons
    """

    if dark_mode:
        plt.style.use('dark_background')

    species_list = get_species_list(A_dict_dNdS)
    print(f"... plotting {len(species_list)} species")

    species_count = len(species_list)
    assert species_count==2
    species_index = {species : i for i, species in enumerate(species_list)}

    cols = species_count
    rows = 1
    fig, axes = plt.subplots(rows, cols, figsize=(16, 10)) # for less than three rows
    fs = 35
    plt.rcParams['text.usetex'] = True

    colors_dict = {
        "grey" : "#3A4040", 
        "bars" : "#7BB7AE", 
        "line" : "#FE4894",
        "pdf" : "#3D7068", 
        "A" : "#F2933A", # uniform_filtered orange
        "X" : "#b82946", # native red
    }

    pairs_sorted = sorted(list(A_dict_dNdS.keys()))
    assert len(pairs_sorted) ==1
    for pair in pairs_sorted:
        ### get pair indices for species pair
        try:
            gen1, spec1, gen2, spec2 =pair.split("_")
        except:
            raise RuntimeError(f"{pair} could not be parsed")
        species1 = f"{gen1}_{spec1}"
        species2 = f"{gen2}_{spec2}"

        species1_lab = species1.replace("_", ". ")
        species2_lab = species2.replace("_", ". ")
        fig.suptitle(f"{species1_lab} vs. {species2_lab}: \nA and X dNdS violin plot and pos. selected genes", fontsize = fs*1.5)
        
        ### plot binary bar charts

        ## exclude all the NaNs because violinplot can't handle them
        data_A_nan = np.array(pos_list_A[pair], dtype=float)
        data_X_nan = np.array(pos_list_X[pair], dtype=float)
        data_A = [dNdS_A for dNdS_A in data_A_nan if not np.isnan(dNdS_A) ]
        data_X = [dNdS_X for dNdS_X in data_X_nan if not np.isnan(dNdS_X) ]
    
        if len(data_A)==0 or len(data_X)==0:
            # axes[col,row].axis('off')
            # axes[col,row].set_title(f'{species1}\n{species2}', fontsize = fs*0.85)
            continue

        n_A = len(data_A)
        n_X = len(data_X)
        mean_A = np.nanmedian(data_A)
        mean_X = np.nanmedian(data_X)

        data_AX = [data_A, data_X]

        site_classes.binary_barplot_single_pair(data_A=data_A, data_X=data_X, row=1, n_A=n_A, n_X=n_X, axes=axes, colors_dict=colors_dict, fs = fs, ylab ="pos. sel. genes")
        print(f"{species1_lab} vs. {species2_lab} --> sig.pos.sel bar plot")

        ### plot violins
        ## exclude all the NaNs because violinplot can't handle them
        data_A_nan = np.array(A_dict_dNdS[pair], dtype=float)
        data_X_nan = np.array(X_dict_dNdS[pair], dtype=float)
        data_A = [dNdS_A for dNdS_A in data_A_nan if not np.isnan(dNdS_A) ]
        data_X = [dNdS_X for dNdS_X in data_X_nan if not np.isnan(dNdS_X) ]
    
        if len(data_A)==0 or len(data_X)==0:
            # axes[col,row].axis('off')
            # axes[col,row].set_title(f'{species1}\n{species2}', fontsize = fs*0.85)
            continue

        n_A = len(data_A)
        n_X = len(data_X)
        mean_A = np.nanmedian(data_A)
        mean_X = np.nanmedian(data_X)

        data_AX = [data_A, data_X]
        
        violinplot_pair_single(data_A_X=data_AX, col=0, n_A=n_A, n_X=n_X, mean_A=mean_A, mean_X=mean_X, axes = axes, colors_dict=colors_dict, fs = fs, xlab = violin_label, ymax = violin_ymax)

        # axes[col,row].set_title(f'{species2}\n{species1}', fontsize = fs*0.85)
        # axes[col,row].set_title(f'{species2}', fontsize = fs)
        print(f"{species1_lab} vs. {species2_lab} --> dNdS violin plot\n")        

    # Adjust layout to prevent overlap (left, bottom, right, top)
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    if dark_mode:
        filename = filename.replace(".png", "_darkmode.png")

    # transparent background
    plt.savefig(filename, dpi = 300, transparent = True)
    # non-transparent background
    filename_tr = filename.replace(".png", "_white_bg.png")
    plt.savefig(filename_tr, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as: {filename} and {filename_tr}")



if __name__ == """__main__""":

    # username = "milena"
    username = "miltr339"
    chromosome = "A"
    data_files = {"A" : ["A_dNdS", "A_LRT"],
                  "X" : ["X_dNdS", "X_LRT"]}
    summary_paths = get_summary_paths(username=username)

    # bruchini
    if True:
        species_excl = ["D_carinulata", "D_sublineata", "T_castaneum", "T_freemani", "C_septempunctata", "C_magnifica"]
        filename =f"/Users/{username}/work/PhD_code/PhD_chapter3/data/fastX_ortholog_ident/fastX_permutation_bruchini.png"
    # coccinella
    elif False:
        species_excl = ["D_carinulata", "D_sublineata", "T_castaneum", "T_freemani", "B_siliquastri", "A_obtectus", "C_maculatus", "C_chinensis"]
        filename =f"/Users/{username}/work/PhD_code/PhD_chapter3/data/fastX_ortholog_ident/fastX_permutation_coccinella.png"
    # tribolium
    elif True:
        species_excl = ["D_carinulata", "D_sublineata", "C_septempunctata", "C_magnifica", "B_siliquastri", "A_obtectus", "C_maculatus", "C_chinensis"]
        filename =f"/Users/{username}/work/PhD_code/PhD_chapter3/data/fastX_ortholog_ident/fastX_permutation_tribolium.png"
    
    print(" * reading dNdS_dict_A ... ")
    dNdS_dict_A = read_filtered_dNdS_summary(summary_paths[data_files["A"][0]], excl_list=species_excl, max_dS=2, max_dNdS=2)
    print(" * reading dNdS_dict_X ... ")
    dNdS_dict_X = read_filtered_dNdS_summary(summary_paths[data_files["X"][0]], excl_list=species_excl, max_dS=2, max_dNdS=2)
    
    species = get_species_list(dNdS_dict_A, exclude_list=species_excl)
    pairs_list = list(dNdS_dict_A.keys())

    bootstraps = { pair : [] for pair in pairs_list}
    median_diffs = {pair : np.NaN for pair in pairs_list}

    
    ### test with 100, takes a bit of time otherwise
    ### actual analysis with 10000
    num_permutations = 10000

    for pair in pairs_list:
        dNdS_A = dNdS_dict_A[pair]
        dNdS_X = dNdS_dict_X[pair]
        median_diffs[pair] = np.nanmedian(dNdS_A) - np.nanmedian(dNdS_X)
        bootstraps[pair] = permutate_dNdS(dNdS_A=dNdS_A, dNdS_X=dNdS_X, num_permut=num_permutations)
        mean_boot = np.mean(bootstraps[pair])
        print(f" *  {pair} --> \t median(dNdS_A)-median(dNdS_X) = {median_diffs[pair]:.3f}, mean bootstrap diff = {mean_boot:.5f}")
    
    if len(pairs_list)>1:
        print(f"plot pairwise matrix ...")
        plot_dNdS_permutations(boot_diff=bootstraps,measure_diff=median_diffs, A_dict=dNdS_dict_A, X_dict=dNdS_dict_X, filename = filename)
    else:
        print(f"plot one pair ...")
        plot_dNdS_permutations_one_pair(boot_diff=bootstraps,measure_diff=median_diffs, A_dict=dNdS_dict_A, X_dict=dNdS_dict_X, filename = filename)